using System.Collections.Generic;

namespace Neuronic.Filters.Chebyshev
{
    public class HighPassChebyshevCoefficients : ChebyshevCoefficients
    {
        public HighPassChebyshevCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
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

            Helpers.HighPassTransform(CutoffFrequency / SamplingFrequency, DigitalProto, AnalogProto);

            DigitalProto.SetLayout(coeffs);

            return 1d;
        }
    }
}

