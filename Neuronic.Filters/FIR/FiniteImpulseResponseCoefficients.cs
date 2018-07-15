using System;
using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    public abstract class FiniteImpulseResponseCoefficients
    {
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

        public TapChain Calculate()
        {
            var coeffs = new List<double>(FilterOrder);
            Calculate(coeffs);
            return new TapChain(coeffs);
        }

        public abstract void Calculate(IList<double> coeffs);
    }
}