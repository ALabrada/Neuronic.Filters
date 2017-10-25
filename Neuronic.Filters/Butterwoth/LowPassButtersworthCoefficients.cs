using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.Filters.Butterwoth
{
    /// <summary>
    /// A low-pass butterworth filter.
    /// </summary>
    public class LowPassButtersworthCoefficients : ButtersworthCoefficients
    {
        private readonly double _freq;

        /// <summary>
        /// Initializes a new instance of <see cref="LowPassButtersworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cut-off frequency.</param>
        public LowPassButtersworthCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
        {
            CutoffFrequency = cutoffFrequency;
            _freq = 2 * Math.Tan(Math.PI * cutoffFrequency / fs);
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <inheritdoc />
        protected override int NumFilters => (FilterOrder + 1) / 2;

        /// <summary>
        ///     Convert analog lowpass prototype poles to lowpass
        /// </summary>
        private static double ConvertToLowPass(double freq, IList<Complex> poles, IList<Complex> zeros)
        {
            var gain = Math.Pow(freq, poles.Count);

            zeros.Clear();
            for (var i = 0; i < poles.Count; i++)
                poles[i] = freq * poles[i];

            return gain;
        }

        /// <inheritdoc />
        protected override double ConvertPoles()
        {
            return ConvertToLowPass(_freq, Poles, Zeros);
        }

        /// <inheritdoc />
        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = preBLTgain * (preBLTgain / gain);
        }
    }
}