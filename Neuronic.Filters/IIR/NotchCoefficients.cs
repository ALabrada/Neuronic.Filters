using System;
using System.Collections.Generic;
using System.Text;

namespace Neuronic.Filters.IIR
{
    public class NotchCoefficients : IBiquadCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="NotchCoefficients"/> class.
        /// </summary>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="cutoffFrequency">The cutoff frequency.</param>
        /// <param name="quality">The quality.</param>
        public NotchCoefficients(double fs, double cutoffFrequency, double quality)
        {
            SamplingFrequency = fs;
            CutoffFrequency = cutoffFrequency;
            Quality = quality;
        }


        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <summary>
        /// Gets the quality.
        /// </summary>
        public double Quality { get; }

        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>
        /// The gain of the filter.
        /// </returns>
        public double Calculate(IList<Biquad> coeffs)
        {
            var w0 = 2 * CutoffFrequency / SamplingFrequency;
            var bw = w0 / Quality;

            bw *= Math.PI;
            w0 *= Math.PI;

            var gb = 1 / Math.Sqrt(2);
            var beta = (Math.Sqrt(1.0 - gb * gb) / gb) * Math.Tan(bw / 2.0);
            var gain = 1.0 / (1.0 + beta);

            coeffs.Add(new Biquad(1, -2.0 * Math.Cos(w0), 1, 1, -2.0 * gain * Math.Cos(w0), 2.0 * gain - 1.0));

            return 1d;
        }

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>(1);
            var gain = Calculate(coeffs);
            return new TransposedDirectFormIIBiquadChain(coeffs, gain);
        }
    }

    public class PeakCoefficients : IBiquadCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="NotchCoefficients"/> class.
        /// </summary>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="peakFrequency">The peak frequency.</param>
        /// <param name="quality">The quality.</param>
        public PeakCoefficients(double fs, double peakFrequency, double quality)
        {
            SamplingFrequency = fs;
            PeakFrequency = peakFrequency;
            Quality = quality;
        }


        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double PeakFrequency { get; }

        /// <summary>
        /// Gets the quality.
        /// </summary>
        public double Quality { get; }

        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>
        /// The gain of the filter.
        /// </returns>
        public double Calculate(IList<Biquad> coeffs)
        {
            var w0 = 2 * PeakFrequency / SamplingFrequency;
            var bw = w0 / Quality;

            bw *= Math.PI;
            w0 *= Math.PI;

            var gb = 1 / Math.Sqrt(2);
            var beta = (gb / Math.Sqrt(1.0 - gb * gb)) * Math.Tan(bw / 2.0);
            var gain = 1.0 / (1.0 + beta);

            coeffs.Add(new Biquad(1, 0, -1, 1, -2.0 * gain * Math.Cos(w0), 2.0 * gain - 1.0));

            return 1.0 - gain;
        }

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>(1);
            var gain = Calculate(coeffs);
            return new TransposedDirectFormIIBiquadChain(coeffs, gain);
        }
    }
}
