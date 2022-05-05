using Neuronic.Filters.IIR;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Neuronic.Filters
{
    internal static class Helpers
    {
        public const double Ln10 = 2.3025850929940456840179914546844;

        public static double Asinh(double x) => Math.Log(x + Math.Sqrt(x * x + 1));

        public static Complex Infinity() => new Complex(double.PositiveInfinity, 0d);

        public static bool IsInfinity(this Complex c) => double.IsInfinity(c.Real) || double.IsInfinity(c.Imaginary);

        public static double MagnitudeSquared(this Complex c) => c.Real * c.Real + c.Imaginary * c.Imaginary;

        public static Complex AddMul(Complex c, double v, Complex c1)
        {
            return new Complex(c.Real + v * c1.Real, c.Imaginary + v * c1.Imaginary);
        }

        public static Complex GetResponse(this IEnumerable<Biquad> stages, double normalizedFrequency)
        {
            double w = 2 * Math.PI * normalizedFrequency;

            var czn1 = Complex.FromPolarCoordinates(1.0, -w);
            var czn2 = Complex.FromPolarCoordinates(1.0, -2 * w);
            var ch = Complex.One;
            var cbot = Complex.One;

            foreach (var stage in stages)
            {
                var cb = Complex.One;
                var ct = new Complex(stage.B0 / stage.A0, 0d);
                ct = AddMul(ct, stage.B1 / stage.A0, czn1);
                ct = AddMul(ct, stage.B2 / stage.A0, czn2);
                cb = AddMul(cb, stage.A1 / stage.A0, czn1);
                cb = AddMul(cb, stage.A2 / stage.A0, czn2);
                ch *= ct;
                cbot *= cb;
            }

            return ch / cbot;
        }

        public static void LowPassTransform(double fc, Layout digital, Layout analog)
        {
            digital.Clear();

            // prewarp
            var f = Math.Tan(Math.PI * fc);

            Complex transform(Complex c)
            {
                if (c.IsInfinity())
                    return new Complex(-1, 0);

                // frequency transform
                c = f * c;

                // bilinear low pass transform
                return (1.0 + c) / (1.0 - c);
            }

            var numPoles = analog.NumberOfPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var pair = analog[i];
                digital.AddPoleZeroConjugatePairs(transform(pair.Poles.First),
                    transform(pair.Zeros.First));
            }

            if ((numPoles & 1) != 0)
            {
                var pair = analog[pairs];
                digital.Add(transform(pair.Poles.First),
                    transform(pair.Zeros.First));
            }

            digital.NormalW = analog.NormalW;
            digital.NormalGain = analog.NormalGain;
        }

        public static void HighPassTransform(double fc, Layout digital, Layout analog)
        {
            digital.Clear();

            // prewarp
            var f = 1.0 / Math.Tan(Math.PI * fc);

            Complex transform(Complex c)
            {
                if (c.IsInfinity())
                    return new Complex(1, 0);

                // frequency transform
                c = f * c;

                // bilinear low pass transform
                return - (1.0 + c) / (1.0 - c);
            }

            var numPoles = analog.NumberOfPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var pair = analog[i];
                digital.AddPoleZeroConjugatePairs(transform(pair.Poles.First),
                    transform(pair.Zeros.First));
            }

            if ((numPoles & 1) != 0)
            {
                var pair = analog[pairs];
                digital.Add(transform(pair.Poles.First),
                    transform(pair.Zeros.First));
            }

            digital.NormalW = Math.PI - analog.NormalW;
            digital.NormalGain = analog.NormalGain;
        }

        public static void BandPassTransform(double fc, double fw, Layout digital, Layout analog)
        {
            digital.Clear();

            var ww = 2 * Math.PI * fw;

            // pre-calcs
            var wc2 = 2 * Math.PI * fc - (ww / 2);
            var wc = wc2 + ww;

            // what is this crap?
            if (wc2 < 1e-8)
                wc2 = 1e-8;
            if (wc > Math.PI - 1e-8)
                wc = Math.PI - 1e-8;

            var a = Math.Cos((wc + wc2) * 0.5) /
                    Math.Cos((wc - wc2) * 0.5);
            var b = 1 / Math.Tan((wc - wc2) * 0.5);
            var a2 = a * a;
            var b2 = b * b;
            var ab = a * b;
            var ab_2 = 2 * ab;

            ComplexPair transform(Complex c)
            {
                if (c.IsInfinity())
                    return new ComplexPair(-1, 1);

                c = (1.0 + c) / (1.0 - c); // bilinear

                var v = Complex.Zero;
                v = AddMul(v, 4 * (b2 * (a2 - 1) + 1), c);
                v += 8 * (b2 * (a2 - 1) - 1);
                v *= c;
                v += 4 * (b2 * (a2 - 1) + 1);
                v = Complex.Sqrt(v);

                var u = -v;
                u = AddMul(u, ab_2, c);
                u += ab_2;

                v = AddMul(v, ab_2, c);
                v += ab_2;

                var d = Complex.Zero;
                d = AddMul(d, 2 * (b - 1), c) + 2 * (1 + b);

                return new ComplexPair(u / d, v / d);
            }

            var numPoles = analog.NumberOfPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var pair = analog[i];
                var p1 = transform(pair.Poles.First);
                var z1 = transform(pair.Zeros.First);

                //
                // Optimize out the calculations for conjugates for Release builds
                //
#if DEBUG
                ComplexPair p2 = transform(pair.Poles.Second);
                ComplexPair z2 = transform(pair.Zeros.Second);
                Debug.Assert(p2.First == Complex.Conjugate(p1.First));
                Debug.Assert(p2.Second == Complex.Conjugate(p1.Second));
#endif

                digital.AddPoleZeroConjugatePairs(p1.First, z1.First);
                digital.AddPoleZeroConjugatePairs(p1.Second, z1.Second);
            }

            if ((numPoles & 1) != 0)
            {
                var poles = transform(analog[pairs].Poles.First);
                var zeros = transform(analog[pairs].Zeros.First);

                digital.Add(poles, zeros);
            }

            double wn = analog.NormalW;
            digital.NormalW = 2 * Math.Atan(Math.Sqrt(Math.Tan((wc + wn) * 0.5) * Math.Tan((wc2 + wn) * 0.5)));
            digital.NormalGain = analog.NormalGain;
        }

        public static void BandStopTransform(double fc, double fw, Layout digital, Layout analog)
        {
            digital.Clear();

            var ww = 2 * Math.PI * fw;

            var wc2 = 2 * Math.PI * fc - (ww / 2);
            var wc = wc2 + ww;

            // this is crap
            if (wc2 < 1e-8)
                wc2 = 1e-8;
            if (wc > Math.PI - 1e-8)
                wc = Math.PI - 1e-8;

            var a = Math.Cos((wc + wc2) * 0.5) /
                    Math.Cos((wc - wc2) * 0.5);
            var b = Math.Tan((wc - wc2) * 0.5);
            var a2 = a * a;
            var b2 = b * b;

            ComplexPair transform(Complex c)
            {
                if (c.IsInfinity())
                    c = -1;
                else
                    c = (1.0 + c) / (1.0 - c); // bilinear

                var u = Complex.Zero;
                u = AddMul(u, 4 * (b2 + a2 - 1), c);
                u += 8 * (b2 - a2 + 1);
                u *= c;
                u += 4 * (a2 + b2 - 1);
                u = Complex.Sqrt(u);

                var v = u * -.5;
                v += a;
                v = AddMul(v, -a, c);

                u *= .5;
                u += a;
                u = AddMul(u, -a, c);

                Complex d = b + 1;
                d = AddMul(d, b - 1, c);

                return new ComplexPair(u / d, v / d);
            }

            var numPoles = analog.NumberOfPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var pair = analog[i];
                var p = transform(pair.Poles.First);
                var z = transform(pair.Zeros.First);

                //
                // Optimize out the calculations for conjugates for Release builds
                //
                if (z.Second == z.First)
                    z.Second = Complex.Conjugate(z.First);

                digital.AddPoleZeroConjugatePairs(p.First, z.First);
                digital.AddPoleZeroConjugatePairs(p.Second, z.Second);
            }

            if ((numPoles & 1) != 0)
            {
                var poles = transform(analog[pairs].Poles.First);
                var zeros = transform(analog[pairs].Zeros.First);

                digital.Add(poles, zeros);
            }

            if (fc < 0.25)
                digital.NormalW = Math.PI;
            else
                digital.NormalW = 0;
            digital.NormalGain = analog.NormalGain;
        }

        public static double SetLayout(this Layout proto, IList<Biquad> stages)
        {
            var numPoles = proto.NumberOfPoles;
            var numStages = (numPoles + 1) / 2;

            for (int i = 0; i < numStages; ++i)
                stages.Insert(0, Biquad.FromPoleZeroPair(proto[i]));

            return proto.NormalGain /
                stages.GetResponse(proto.NormalW / (2 * Math.PI)).Magnitude;
        }

        public static void Scale(this IList<Biquad> stages, double scale)
        {
            stages[0] *= scale;
        }

        //public static double ToZPK(this IList<double> b, IList<Complex> z)
        //{
        //    var k = b[0];

        //    foreach (var root in b.FindRoots())
        //        z.Add(root);

        //    return k;
        //}

        //public static IEnumerable<Complex> FindRoots(this IList<double> polynomial)
        //{
        //    var order = polynomial.Count - 1;
        //    var matrix = new double[order, order];
        //    var den = polynomial[0];
        //    for (int i = 0; i < order - 1; i++)
        //    {
        //        matrix[i, i + 1] = 1;
        //        matrix[i, 0] = -polynomial[i + 1] / den;
        //    }
        //    matrix[order - 1, 0] = -polynomial[order] / den;

        //    var d = new double[order];
        //    var e = new double[order];

        //    for (int i = 0; i < order; i++)
        //        d[i] = matrix[order - 1, i];

        //    Evd(matrix, d, e, order);

        //    for (int i = 0; i < order; i++)
        //        yield return new Complex(d[i], e[i]);
        //}

        //// Adapted from MathNet.Numerics
        //internal static void Evd(double[,] matrix, double[] d, double[] e, int order)
        //{
        //    unsafe
        //    {
        //        fixed (double* aPtr = matrix, dPtr = d, ePtr = e)
        //        {
        //            SymmetricTridiagonalize(aPtr, dPtr, ePtr, order);
        //            SymmetricDiagonalize(aPtr, dPtr, ePtr, order);
        //        }
        //    }

        //}

        //// Adapted from MathNet.Numerics
        //internal static unsafe void SymmetricTridiagonalize(double* a, double* d, double* e, int order)
        //{
        //    // Householder reduction to tridiagonal form.
        //    for (var i = order - 1; i > 0; i--)
        //    {
        //        // Scale to avoid under/overflow.
        //        var scale = 0.0;
        //        var h = 0.0;

        //        for (var k = 0; k < i; k++)
        //        {
        //            scale = scale + Math.Abs(d[k]);
        //        }

        //        if (scale == 0.0)
        //        {
        //            e[i] = d[i - 1];
        //            for (var j = 0; j < i; j++)
        //            {
        //                d[j] = a[(j * order) + i - 1];
        //                a[(j * order) + i] = 0.0;
        //                a[(i * order) + j] = 0.0;
        //            }
        //        }
        //        else
        //        {
        //            // Generate Householder vector.
        //            for (var k = 0; k < i; k++)
        //            {
        //                d[k] /= scale;
        //                h += d[k] * d[k];
        //            }

        //            var f = d[i - 1];
        //            var g = Math.Sqrt(h);
        //            if (f > 0)
        //            {
        //                g = -g;
        //            }

        //            e[i] = scale * g;
        //            h = h - (f * g);
        //            d[i - 1] = f - g;

        //            for (var j = 0; j < i; j++)
        //            {
        //                e[j] = 0.0;
        //            }

        //            // Apply similarity transformation to remaining columns.
        //            for (var j = 0; j < i; j++)
        //            {
        //                f = d[j];
        //                a[(i * order) + j] = f;
        //                g = e[j] + (a[(j * order) + j] * f);

        //                for (var k = j + 1; k <= i - 1; k++)
        //                {
        //                    g += a[(j * order) + k] * d[k];
        //                    e[k] += a[(j * order) + k] * f;
        //                }

        //                e[j] = g;
        //            }

        //            f = 0.0;

        //            for (var j = 0; j < i; j++)
        //            {
        //                e[j] /= h;
        //                f += e[j] * d[j];
        //            }

        //            var hh = f / (h + h);

        //            for (var j = 0; j < i; j++)
        //            {
        //                e[j] -= hh * d[j];
        //            }

        //            for (var j = 0; j < i; j++)
        //            {
        //                f = d[j];
        //                g = e[j];

        //                for (var k = j; k <= i - 1; k++)
        //                {
        //                    a[(j * order) + k] -= (f * e[k]) + (g * d[k]);
        //                }

        //                d[j] = a[(j * order) + i - 1];
        //                a[(j * order) + i] = 0.0;
        //            }
        //        }

        //        d[i] = h;
        //    }

        //    // Accumulate transformations.
        //    for (var i = 0; i < order - 1; i++)
        //    {
        //        a[(i * order) + order - 1] = a[(i * order) + i];
        //        a[(i * order) + i] = 1.0;
        //        var h = d[i + 1];
        //        if (h != 0.0)
        //        {
        //            for (var k = 0; k <= i; k++)
        //            {
        //                d[k] = a[((i + 1) * order) + k] / h;
        //            }

        //            for (var j = 0; j <= i; j++)
        //            {
        //                var g = 0.0;
        //                for (var k = 0; k <= i; k++)
        //                {
        //                    g += a[((i + 1) * order) + k] * a[(j * order) + k];
        //                }

        //                for (var k = 0; k <= i; k++)
        //                {
        //                    a[(j * order) + k] -= g * d[k];
        //                }
        //            }
        //        }

        //        for (var k = 0; k <= i; k++)
        //        {
        //            a[((i + 1) * order) + k] = 0.0;
        //        }
        //    }

        //    for (var j = 0; j < order; j++)
        //    {
        //        d[j] = a[(j * order) + order - 1];
        //        a[(j * order) + order - 1] = 0.0;
        //    }

        //    a[(order * order) - 1] = 1.0;
        //    e[0] = 0.0;
        //}

        //// Adapted from MathNet.Numerics
        //internal static unsafe void SymmetricDiagonalize(double* a, double* d, double* e, int order)
        //{
        //    const int maxiter = 1000;

        //    double Hypotenuse(double x, double y)
        //    {
        //        if (Math.Abs(x) > Math.Abs(y))
        //        {
        //            double r = y / x;
        //            return Math.Abs(x) * Math.Sqrt(1 + (r * r));
        //        }

        //        if (y != 0.0)
        //        {
        //            double r = x / y;
        //            return Math.Abs(y) * Math.Sqrt(1 + (r * r));
        //        }

        //        return 0d;
        //    }

        //    for (var i = 1; i < order; i++)
        //    {
        //        e[i - 1] = e[i];
        //    }

        //    e[order - 1] = 0.0;

        //    var f = 0.0;
        //    var tst1 = 0.0;
        //    var eps = double.Epsilon;
        //    for (var l = 0; l < order; l++)
        //    {
        //        // Find small subdiagonal element
        //        tst1 = Math.Max(tst1, Math.Abs(d[l]) + Math.Abs(e[l]));
        //        var m = l;
        //        while (m < order)
        //        {
        //            if (Math.Abs(e[m]) <= eps * tst1)
        //            {
        //                break;
        //            }

        //            m++;
        //        }

        //        // If m == l, d[l] is an eigenvalue,
        //        // otherwise, iterate.
        //        if (m > l)
        //        {
        //            var iter = 0;
        //            do
        //            {
        //                iter = iter + 1; // (Could check iteration count here.)

        //                // Compute implicit shift
        //                var g = d[l];
        //                var p = (d[l + 1] - g) / (2.0 * e[l]);
        //                var r = Hypotenuse(p, 1.0);
        //                if (p < 0)
        //                {
        //                    r = -r;
        //                }

        //                d[l] = e[l] / (p + r);
        //                d[l + 1] = e[l] * (p + r);

        //                var dl1 = d[l + 1];
        //                var h = g - d[l];
        //                for (var i = l + 2; i < order; i++)
        //                {
        //                    d[i] -= h;
        //                }

        //                f = f + h;

        //                // Implicit QL transformation.
        //                p = d[m];
        //                var c = 1.0;
        //                var c2 = c;
        //                var c3 = c;
        //                var el1 = e[l + 1];
        //                var s = 0.0;
        //                var s2 = 0.0;
        //                for (var i = m - 1; i >= l; i--)
        //                {
        //                    c3 = c2;
        //                    c2 = c;
        //                    s2 = s;
        //                    g = c * e[i];
        //                    h = c * p;
        //                    r = Hypotenuse(p, e[i]);
        //                    e[i + 1] = s * r;
        //                    s = e[i] / r;
        //                    c = p / r;
        //                    p = (c * d[i]) - (s * g);
        //                    d[i + 1] = h + (s * ((c * g) + (s * d[i])));

        //                    // Accumulate transformation.
        //                    for (var k = 0; k < order; k++)
        //                    {
        //                        h = a[((i + 1) * order) + k];
        //                        a[((i + 1) * order) + k] = (s * a[(i * order) + k]) + (c * h);
        //                        a[(i * order) + k] = (c * a[(i * order) + k]) - (s * h);
        //                    }
        //                }

        //                p = (-s) * s2 * c3 * el1 * e[l] / dl1;
        //                e[l] = s * p;
        //                d[l] = c * p;

        //                // Check for convergence. If too many iterations have been performed,
        //                // throw exception that Convergence Failed
        //                if (iter >= maxiter)
        //                {
        //                    throw new Exception();
        //                }
        //            } while (Math.Abs(e[l]) > eps * tst1);
        //        }

        //        d[l] = d[l] + f;
        //        e[l] = 0.0;
        //    }

        //    // Sort eigenvalues and corresponding vectors.
        //    for (var i = 0; i < order - 1; i++)
        //    {
        //        var k = i;
        //        var p = d[i];
        //        for (var j = i + 1; j < order; j++)
        //        {
        //            if (d[j] < p)
        //            {
        //                k = j;
        //                p = d[j];
        //            }
        //        }

        //        if (k != i)
        //        {
        //            d[k] = d[i];
        //            d[i] = p;
        //            for (var j = 0; j < order; j++)
        //            {
        //                p = a[(i * order) + j];
        //                a[(i * order) + j] = a[(k * order) + j];
        //                a[(k * order) + j] = p;
        //            }
        //        }
        //    }
        //}
    }
}
