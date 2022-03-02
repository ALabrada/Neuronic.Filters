using Neuronic.Filters.IIR;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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
    }
}
