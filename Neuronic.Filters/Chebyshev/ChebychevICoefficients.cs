using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace Neuronic.Filters.Chebyshev
{
    public abstract class ChebyshevICoefficients : IBiquadCoefficients
    {
        protected ChebyshevICoefficients(int filterOrder, double fs, double rippleDb)
        {
            FilterOrder = filterOrder;
            SamplingFrequency = fs;
            RippleDb = rippleDb;

            AnalogProto = new Layout(FilterOrder);
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }

        public double RippleDb { get; }

        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        internal Layout AnalogProto { get; }

        protected void AnalogDesign()
        {
            var numPoles = FilterOrder;
            var rippleDb = RippleDb;

            AnalogProto.Clear();

            var eps = Math.Sqrt(1.0 / Math.Exp(-rippleDb * 0.1 * Helpers.Ln10) - 1);
            var v0 = Helpers.Asinh(1 / eps) / numPoles;
            var sinh_v0 = -Math.Sinh(v0);
            var cosh_v0 = Math.Cosh(v0);

            var n2 = 2 * numPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var k = 2 * i + 1 - numPoles;
                var a = sinh_v0 * Math.Cos(k * Math.PI / n2);
                var b = cosh_v0 * Math.Sin(k * Math.PI / n2);

                //addPoleZero (complex_t (a, b), infinity());
                //addPoleZero (complex_t (a, -b), infinity());
                AnalogProto.AddPoleZeroConjugatePairs(new Complex(a, b), Helpers.Infinity());
            }

            if ((numPoles & 1) != 0)
            {
                AnalogProto.Add(new Complex(sinh_v0, 0), Helpers.Infinity());
                AnalogProto.NormalW = 0;
                AnalogProto.NormalGain = 1;
            }
            else
            {
                AnalogProto.NormalW = 0;
                AnalogProto.NormalGain = Math.Pow(10, -rippleDb / 20.0);
            }
        }

        public abstract double Calculate(IList<Biquad> coeffs);

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public virtual BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>((FilterOrder + 1) / 2);
            var gain = Calculate(coeffs);
            return new TransposedDirectFormIIBiquadChain(coeffs, gain);
        }
    }
}

