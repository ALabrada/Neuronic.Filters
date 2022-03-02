using Neuronic.Filters.IIR;
using System.Numerics;

namespace Neuronic.Filters
{
    /// <summary>
    /// Structure that holds the coefficients of a biquad filter (see <see ref="https://en.wikipedia.org/wiki/Digital_biquad_filter"/>).
    /// </summary>
    public struct Biquad
    {
        private readonly double _a1, _a2;
        private double _b0, _b1, _b2;

        /// <summary>
        /// Initializes a new instance of <see cref="Biquad"/>.
        /// </summary>
        /// <param name="b0">The b<sub>0</sub> coefficient.</param>
        /// <param name="b1">The b<sub>1</sub> coefficient.</param>
        /// <param name="b2">The b<sub>2</sub> coefficient.</param>
        /// <param name="a0">The a<sub>0</sub> coefficient.</param>
        /// <param name="a1">The a<sub>1</sub> coefficient.</param>
        /// <param name="a2">The a<sub>2</sub> coefficient.</param>
        public Biquad(double b0, double b1, double b2, double a0, double a1, double a2)
        {
            _b0 = b0 / a0;
            _b1 = b1 / a0;
            _b2 = b2 / a0;
            _a1 = a1 / a0;
            _a2 = a2 / a0;
        }
        /// <summary>
        /// The b<sub>0</sub> coefficient.
        /// </summary>
        public double B0 => _b0;
        /// <summary>
        /// The b<sub>1</sub> coefficient.
        /// </summary>
        public double B1 => _b1;
        /// <summary>
        /// The b<sub>2</sub> coefficient.
        /// </summary>
        public double B2 => _b2;
        /// <summary>
        /// The a<sub>0</sub> coefficient.
        /// </summary>
        public double A0 => 1d;
        /// <summary>
        /// The a<sub>1</sub> coefficient.
        /// </summary>
        public double A1 => _a1;
        /// <summary>
        /// The a<sub>2</sub> coefficient.
        /// </summary>
        public double A2 => _a2;

        public static readonly Biquad Identity = new Biquad(1, 0, 0, 1, 0, 0);

        internal static Biquad FromPoleZeroPair(IIR.PoleZeroPair pair)
        {
            if (pair.IsSignlePole)
                return FromPole(pair.Poles.First, pair.Zeros.First);
            else
                return FromPoles(
                    pair.Poles.First, pair.Zeros.First,
                    pair.Poles.Second, pair.Zeros.Second);
        }

        public static Biquad FromPole(Complex pole, Complex zero)
        {
            return new Biquad(-zero.Real, 1, 0, 1, -pole.Real, 0);
        }

        public static Biquad FromPoles(
            Complex pole1, Complex zero1,
            Complex pole2, Complex zero2)
        {
            const double a0 = 1;
            double a1;
            double a2;

            if (pole1.Imaginary != 0)
            {
                a1 = -2 * pole1.Real;
                a2 = pole1.MagnitudeSquared();
            }
            else
            {
                a1 = -(pole1.Real + pole2.Real);
                a2 = pole1.Real * pole2.Real;
            }

            const double b0 = 1;
            double b1;
            double b2;

            if (zero1.Imaginary != 0)
            {
                b1 = -2 * zero1.Real;
                b2 = zero1.MagnitudeSquared();
            }
            else
            {
                b1 = -(zero1.Real + zero2.Real);
                b2 = zero1.Real * zero2.Real;
            }

            return new Biquad(b0, b1, b2, a0, a1, a2);
        }

        public static Biquad operator *(Biquad b, double scale)
        {
            b._b0 *= scale;
            b._b1 *= scale;
            b._b2 *= scale;
            return b;
        }

        /// <summary>Returns the fully qualified type name of this instance.</summary>
        /// <returns>A <see cref="T:System.String" /> containing a fully qualified type name.</returns>
        /// <filterpriority>2</filterpriority>
        public override string ToString()
        {
            return $"B: [{B0:F4}, {B1:F4}, {B2:F4}], A: [{A0:F4}, {A1:F4}, {A2:F4}]";
        }
    }
}