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

        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            Helpers.LowPassTransform(CutoffFrequency / SamplingFrequency, DigitalProto, AnalogProto);

            DigitalProto.SetLayout(coeffs);

            return 1d;
        }
    }
}