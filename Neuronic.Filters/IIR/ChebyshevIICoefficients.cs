using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// Base class for the Inverse Chebyshev filter designers.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IBiquadCoefficients" />
    public abstract class ChebyshevIICoefficients : IBiquadCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="ChebyshevIICoefficients"/> class.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="stopBandDb">The stop band ripple in Decibels.</param>
        protected ChebyshevIICoefficients(int filterOrder, double fs, double stopBandDb)
        {
            FilterOrder = filterOrder;
            SamplingFrequency = fs;
            StopBandDb = stopBandDb;

            AnalogProto = new Layout(FilterOrder);
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }

        /// <summary>
        /// Gets the maximum stopband ripple in Decibels.
        /// </summary>
        public double StopBandDb { get; }

        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        internal Layout AnalogProto { get; }


        /// <summary>
        /// Designs the analog low pass filter prototype.
        /// </summary>
        protected void AnalogDesign()
        {
            var numPoles = FilterOrder;
            var stopBandDb = StopBandDb;

            AnalogProto.Clear();

            var eps = Math.Sqrt(1.0 / (Math.Exp(stopBandDb * 0.1 * Helpers.Ln10) - 1));
            var v0 = Helpers.Asinh(1 / eps) / numPoles;
            var sinh_v0 = -Math.Sinh(v0);
            var cosh_v0 = Math.Cosh(v0);
            var fn = Math.PI / (2 * numPoles);

            var k = 1;
            for (int i = numPoles / 2; --i >= 0; k += 2)
            {
                var a = sinh_v0 * Math.Cos((k - numPoles) * fn);
                var b = cosh_v0 * Math.Sin((k - numPoles) * fn);
                var d2 = a * a + b * b;
                var im = 1 / Math.Cos(k * fn);

                AnalogProto.AddPoleZeroConjugatePairs(
                    new Complex(a / d2, b / d2),
                    new Complex(0, im));
            }

            if ((numPoles & 1) != 0)
            {
                AnalogProto.Add(1 / sinh_v0, Helpers.Infinity());
            }
        }

        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>
        /// The gain of the filter.
        /// </returns>
        public abstract double Calculate(IList<Biquad> coeffs);

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public virtual BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>((FilterOrder + 1) / 2);
            var gain = Calculate(coeffs);
            coeffs.Scale(gain);
            return new TransposedDirectFormIIBiquadChain(coeffs, 1);
        }
    }
}

