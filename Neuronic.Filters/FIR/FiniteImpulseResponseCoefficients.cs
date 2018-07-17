using System;
using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    /// <summary>
    /// Abstraction of a Finite Impulse Response (FIR) filter designer.
    /// </summary>
    public abstract class FiniteImpulseResponseCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="FiniteImpulseResponseCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order..</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// n
        /// or
        /// fs
        /// </exception>
        public FiniteImpulseResponseCoefficients(int n, double fs)
        {
            if (n <= 0)
                throw new ArgumentOutOfRangeException(nameof(n));
            if (fs <= 0)
                throw new ArgumentOutOfRangeException(nameof(fs));

            FilterOrder = n;
            SamplingFrequency = fs;
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }

        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        /// <summary>
        /// Designs a FIR filter.
        /// </summary>
        /// <returns>The FIR filter.</returns>
        public TapChain Calculate()
        {
            var coeffs = new List<double>(FilterOrder);
            Calculate(coeffs);
            return new TapChain(coeffs);
        }

        /// <summary>
        /// Calculates coefficients (taps) for a FIR filter.
        /// </summary>
        /// <param name="coeffs">The collection that will hold the coefficients.</param>
        public abstract void Calculate(IList<double> coeffs);
    }
}