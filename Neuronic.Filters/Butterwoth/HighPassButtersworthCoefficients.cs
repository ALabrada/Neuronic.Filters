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

        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            Helpers.HighPassTransform(CutoffFrequency / SamplingFrequency, DigitalProto, AnalogProto);

            DigitalProto.SetLayout(coeffs);

            return 1d;
        }
    }
}