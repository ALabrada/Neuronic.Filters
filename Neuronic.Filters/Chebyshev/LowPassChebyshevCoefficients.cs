using System.Collections.Generic;

namespace Neuronic.Filters.Chebyshev
{
    public class LowPassChebyshevICoefficients : ChebyshevICoefficients
    {
        public LowPassChebyshevICoefficients(int filterOrder, double fs, double cutoffFrequency, double rippleDb) : base(filterOrder, fs, rippleDb)
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

    public class LowPassChebyshevIICoefficients : ChebyshevIICoefficients
    {
        public LowPassChebyshevIICoefficients(int filterOrder, double fs, double cutoffFrequency, double stopBandDb) : base(filterOrder, fs, stopBandDb)
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

