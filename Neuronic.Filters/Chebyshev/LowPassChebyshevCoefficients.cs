using System.Collections.Generic;

namespace Neuronic.Filters.Chebyshev
{
    public class LowPassChebyshevCoefficients : ChebyshevCoefficients
    {
        public LowPassChebyshevCoefficients(int filterOrder, double fs, double cutoffFrequency) : base(filterOrder, fs)
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

            Helpers.LowPassTransform(CutoffFrequency / SamplingFrequency, DigitalProto, AnalogProto);

            DigitalProto.SetLayout(coeffs);

            return 1d;
        }
    }
}

