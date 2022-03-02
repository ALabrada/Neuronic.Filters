using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// A band-pass butterworth filter.
    /// </summary>
    public class BandPassButterworthCoefficients
        : ButterworthCoefficients

    {
        /// <summary>
        /// Initializes a new instance of <see cref="BandPassButterworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The minor cut-off frequency.</param>
        /// <param name="f2">The major cut-off frequency.</param>
        public BandPassButterworthCoefficients(int filterOrder, double fs, double f1, double f2) : base(filterOrder, fs)
        {
            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }

        /// <inheritdoc />
        protected override int NumFilters => FilterOrder;

        /// <inheritdoc />
        protected override double ConvertPoles()
        {
            var f1 = 2 * Math.Tan(Math.PI * FirstCutoffFrequency / SamplingFrequency);
            var f2 = 2 * Math.Tan(Math.PI * SecondCutoffFrequency / SamplingFrequency);
            return ConvertBandPass(f1, f2, Poles, Zeros);
        }

        /// <inheritdoc />
        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = preBLTgain * (preBLTgain / gain);
        }
    }
}