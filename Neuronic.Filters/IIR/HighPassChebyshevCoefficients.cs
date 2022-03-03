using System.Collections.Generic;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// A high-pass Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevICoefficients" />
    public class HighPassChebyshevICoefficients : ChebyshevICoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="HighPassChebyshevICoefficients"/> class.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cutoff frequency.</param>
        /// <param name="rippleDb">The maximum pass-band ripple.</param>
        public HighPassChebyshevICoefficients(int filterOrder, double fs, double cutoffFrequency, double rippleDb) : base(filterOrder, fs, rippleDb)
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

            Helpers.HighPassTransform(CutoffFrequency / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }

    /// <summary>
    /// A high-pass Inverse Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevIICoefficients" />
    public class HighPassChebyshevIICoefficients : ChebyshevIICoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="HighPassChebyshevIICoefficients"/> class.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cutoff frequency.</param>
        /// <param name="stopBandDb">The maximum stop-band ripple.</param>
        public HighPassChebyshevIICoefficients(int filterOrder, double fs, double cutoffFrequency, double stopBandDb) : base(filterOrder, fs, stopBandDb)
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

            Helpers.HighPassTransform(CutoffFrequency / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }
}

