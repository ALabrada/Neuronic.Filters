﻿using System;
using System.Collections.Generic;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// A band-pass Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevICoefficients" />
    public class BandPassChebyshevICoefficients : ChebyshevICoefficients
    {
        /// <summary>
        /// Initializes a new instance of <see cref="BandPassChebyshevICoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The minor cut-off frequency.</param>
        /// <param name="f2">The major cut-off frequency.</param>
        /// <param name="rippleDb">The maximum passband ripple.</param>
        public BandPassChebyshevICoefficients(int filterOrder, double fs, double f1, double f2, double rippleDb) : base(filterOrder, fs, rippleDb)
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

        /// <inheritdoc/>
        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            var digitalProto = new Layout(AnalogProto.Count * 2);

            var center = (FirstCutoffFrequency + SecondCutoffFrequency) / 2;
            var width = Math.Abs(FirstCutoffFrequency - center);
            Helpers.BandPassTransform(center / SamplingFrequency, width / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }

    /// <summary>
    /// A band-pass Inverse Chebyshev filter designer.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.IIR.ChebyshevIICoefficients" />
    public class BandPassChebyshevIICoefficients : ChebyshevIICoefficients
    {
        /// <summary>
        /// Initializes a new instance of <see cref="BandPassChebyshevIICoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The minor cut-off frequency.</param>
        /// <param name="f2">The major cut-off frequency.</param>
        /// <param name="stopBandDb">The maximum stop-band ripple.</param>
        public BandPassChebyshevIICoefficients(int filterOrder, double fs, double f1, double f2, double stopBandDb) : base(filterOrder, fs, stopBandDb)
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

        /// <inheritdoc/>
        public override double Calculate(IList<Biquad> coeffs)
        {
            AnalogDesign();

            var digitalProto = new Layout(AnalogProto.Count * 2);

            var center = (FirstCutoffFrequency + SecondCutoffFrequency) / 2;
            var width = Math.Abs(FirstCutoffFrequency - center);
            Helpers.BandPassTransform(center / SamplingFrequency, width / SamplingFrequency, digitalProto, AnalogProto);

            return digitalProto.SetLayout(coeffs);
        }
    }
}

