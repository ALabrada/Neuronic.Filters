using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.Butterwoth
{
    /// <summary>
    /// A band-pass butterworth filter.
    /// </summary>
    public class BandPassButtersworthCoefficients : ButtersworthCoefficients
    {
        private readonly double _f1, _f2;

        /// <summary>
        /// Initializes a new instance of <see cref="BandPassButtersworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The minor cut-off frequency.</param>
        /// <param name="f2">The major cut-off frequency.</param>
        public BandPassButtersworthCoefficients(int filterOrder, double fs, double f1, double f2) : base(filterOrder, fs)
        {
            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;
            _f1 = 2 * Math.Tan(Math.PI * f1 / fs);
            _f2 = 2 * Math.Tan(Math.PI * f2 / fs);
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }

        /// <inheritdoc />
        protected override int NumFilters => FilterOrder;

        /// <inheritdoc />
        protected override double ConvertPoles()
        {
            return ConvertBandPass(_f1, _f2, Poles, Zeros);
        }

        /// <inheritdoc />
        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = preBLTgain * (preBLTgain / gain);
        }

        private static double ConvertBandPass(double f1, double f2, IList<Complex> poles,
            IList<Complex> zeros)
        {
            var bw = f2 - f1;
            var wc = Math.Sqrt(f1 * f2);

            // Calculate bandpass gain
            var gain = Math.Pow(bw, poles.Count - zeros.Count);

            // Convert LP poles to BP: these two sets of for-loops result in an ordered
            // list of poles and their complex conjugates
            var tempPoles = new List<Complex>(2 * poles.Count);
            // First set of poles + conjugates
            tempPoles.AddRange(from pole in poles
                where !pole.Equals(Complex.Zero)
                let firstTerm = 0.5 * pole * bw
                let secondTerm = 0.5 * Complex.Sqrt(bw * bw * (pole * pole) - 4 * wc * wc)
                select firstTerm + secondTerm);
            // Second set of poles + conjugates
            tempPoles.AddRange(from pole in poles
                where !pole.Equals(Complex.Zero)
                let firstTerm = 0.5 * pole * bw
                let secondTerm = 0.5 * Complex.Sqrt(bw * bw * (pole * pole) - 4 * wc * wc)
                select firstTerm - secondTerm);
            // Init zeros, no non-zero values to convert
            zeros.Clear();
            for (int i = 0; i < poles.Count; i++)
                zeros.Add(Complex.Zero);
            // Copy converted poles to output array
            poles.Clear();
            foreach (var pole in tempPoles)
                poles.Add(pole);

            return gain;
        }
    }
}