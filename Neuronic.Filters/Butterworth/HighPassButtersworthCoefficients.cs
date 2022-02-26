using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.Butterworth
{
    /// <summary>
    /// A high-pass butterworth filter.
    /// </summary>
    public class HighPassButterworthCoefficients : ButterworthCoefficients
    {
        /// <summary>
        /// Initializes a new instance of <see cref="HighPassButterworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cut-off frequency.</param>
        public HighPassButterworthCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
        {
            CutoffFrequency = cutoffFrequency;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <inheritdoc />
        protected override int NumFilters => (FilterOrder + 1) / 2;

        /// <inheritdoc />
        protected override double ConvertPoles()
        {
            var freq = 2 * Math.Tan(Math.PI * CutoffFrequency / SamplingFrequency);
            return ConvertToHighPass(freq, Poles, Zeros);
        }

        /// <inheritdoc />
        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = 1d / ba[0];
        }
    }
}