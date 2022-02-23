using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.Butterworth
{
    /// <summary>
    /// Base class for the butterwoth filters.
    /// </summary>
    public abstract class ButterworthCoefficients
        : IBiquadCoefficients
    {
        /// <summary>
        /// Creates a new instance of <see cref="ButterworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        protected ButterworthCoefficients(int filterOrder, double fs)
        {
            FilterOrder = filterOrder;
            SamplingFrequency = fs;

            AnalogProto = new Layout(FilterOrder);
            DigitalProto = new Layout(FilterOrder);
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }
        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        internal Layout AnalogProto { get; }

        internal Layout DigitalProto { get; }

        protected void AnalogDesign()
        {
            var numPoles = FilterOrder;

            AnalogProto.Clear();

            var n2 = 2 * numPoles;
            var pairs = numPoles / 2;
            for (int i = 0; i < pairs; ++i)
            {
                var c = Complex.FromPolarCoordinates(1.0, Math.PI / 2.0 + (2 * i + 1) * Math.PI / n2);
                AnalogProto.AddPoleZeroConjugatePairs(c, Helpers.Infinity());
            }

            if ((numPoles & 1) != 0)
                AnalogProto.Add(-1, Helpers.Infinity());
        }

        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>The gain of the filter.</returns>
        public abstract double Calculate(IList<Biquad> coeffs);

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>((FilterOrder + 1) / 2);
            Calculate(coeffs);
            return new DirectFormIBiquadChain(coeffs);
        }
    }
}