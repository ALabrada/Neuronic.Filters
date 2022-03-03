using System.Collections.Generic;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// A low-pass Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevICoefficients" />
    public class LowPassChebyshevICoefficients : ChebyshevICoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="LowPassChebyshevICoefficients"/> class.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cutoff frequency.</param>
        /// <param name="rippleDb">The maximum pass-band ripple.</param>
        public LowPassChebyshevICoefficients(int filterOrder, double fs, double cutoffFrequency, double rippleDb) : base(filterOrder, fs, rippleDb)
        {
            CutoffFrequency = cutoffFrequency;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <inheritdoc/>
        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            var digitalProto = new Layout(AnalogProto.Count);

            Helpers.LowPassTransform(CutoffFrequency / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }

    /// <summary>
    /// A low-pass Inverse Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevIICoefficients" />
    public class LowPassChebyshevIICoefficients : ChebyshevIICoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="LowPassChebyshevIICoefficients"/> class.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cutoff frequency.</param>
        /// <param name="stopBandDb">The maximum stop-band ripple.</param>
        public LowPassChebyshevIICoefficients(int filterOrder, double fs, double cutoffFrequency, double stopBandDb) : base(filterOrder, fs, stopBandDb)
        {
            CutoffFrequency = cutoffFrequency;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <inheritdoc/>
        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            var digitalProto = new Layout(AnalogProto.Count);

            Helpers.LowPassTransform(CutoffFrequency / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }
}

