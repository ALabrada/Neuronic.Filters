using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.Filters.Butterworth
{
    /// <summary>
    /// A low-pass butterworth filter.
    /// </summary>
    public class LowPassButterworthCoefficients
        : ButterworthCoefficients
    {
        /// <summary>
        /// Initializes a new instance of <see cref="LowPassButterworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cut-off frequency.</param>
        public LowPassButterworthCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
        {
            CutoffFrequency = cutoffFrequency;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            var digitalProto = new Layout(AnalogProto.Count);

            Helpers.LowPassTransform(CutoffFrequency / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }
}