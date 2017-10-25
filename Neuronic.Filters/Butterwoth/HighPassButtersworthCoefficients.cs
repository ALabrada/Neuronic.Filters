using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.Butterwoth
{
    /// <summary>
    /// A high-pass butterworth filter.
    /// </summary>
    public class HighPassButtersworthCoefficients : ButtersworthCoefficients
    {
        private readonly double _freq;

        /// <summary>
        /// Initializes a new instance of <see cref="HighPassButtersworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cut-off frequency.</param>
        public HighPassButtersworthCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
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
        /// Convert lowpass poles & zeros to highpass with Wc = f2, use:  hp_S = Wc / lp_S;
        /// </summary>
        /// <param name="freq">The freq.</param>
        /// <param name="poles">The poles.</param>
        /// <param name="zeros">The zeros.</param>
        /// <returns></returns>
        private static double ConvertToHighPass(double freq, IList<Complex> poles, IList<Complex> zeros)
        {
            // Calculate gain
            var prodz = zeros.Aggregate(new Complex(1, 0), (current, zero) => current * -zero);
            var prodp = poles.Aggregate(new Complex(1, 0), (current, pole) => current * -pole);
            var gain = prodz.Real / prodp.Real;

            // Convert LP poles to HP
            for (int i = 0; i < poles.Count; i++)
                if (!poles[i].Equals(Complex.Zero))
                    poles[i] = freq / poles[i];

            // Init with zeros, no non-zero values to convert
            zeros.Clear();
            for (int i = 0; i < poles.Count; i++)
                zeros.Add(Complex.Zero);

            return gain;
        }

        /// <inheritdoc />
        protected override double ConvertPoles()
        {
            return ConvertToHighPass(_freq, Poles, Zeros);
        }

        /// <inheritdoc />
        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = 1d / ba[0];
        }
    }
}