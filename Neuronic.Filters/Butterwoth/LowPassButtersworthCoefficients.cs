using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.Filters.Butterwoth
{
    public class LowPassButtersworthCoefficients : ButtersworthCoefficients
    {
        private readonly double _freq;

        public LowPassButtersworthCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
        {
            CutoffFrequency = cutoffFrequency;
            _freq = 2 * Math.Tan(Math.PI * cutoffFrequency / fs);
        }

        public double CutoffFrequency { get; }

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

        protected override double ConvertPoles()
        {
            return ConvertToLowPass(_freq, Poles, Zeros);
        }

        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = preBLTgain * (preBLTgain / gain);
        }
    }
}