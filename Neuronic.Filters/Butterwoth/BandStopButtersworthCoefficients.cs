using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.Butterwoth
{
    public class BandStopButtersworthCoefficients : ButtersworthCoefficients
    {
        private readonly double _f1, _f2;

        public BandStopButtersworthCoefficients(int filterOrder, double fs, double f1, double f2) : base(filterOrder, fs)
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

        public double FirstCutoffFrequency { get; }
        public double SecondCutoffFrequency { get; }

        protected override int NumFilters => FilterOrder;

        protected override void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
            ba[0] = 1d / ba[0];
        }

        protected override double ConvertPoles()
        {
            return ConvertBandStop(_f1, _f2, Poles, Zeros);
        }

        private static double ConvertBandStop(double f1, double f2, IList<Complex> poles,
            IList<Complex> zeros)
        {
            var bw = f2 - f1;
            var wc = Math.Sqrt(f1 * f2);

            // Compute gain
            var prodz = zeros.Aggregate(new Complex(1, 0), (current, zero) => current * -zero);
            var prodp = poles.Aggregate(new Complex(1, 0), (current, pole) => current * -pole);
            var gain = prodz.Real / prodp.Real;

            // Convert LP zeros to band stop
            var ztmp = new List<Complex>(2 * poles.Count);
            for (int i = 0; i < poles.Count; i++)
            {
                ztmp.Add(new Complex(0, wc));
                ztmp.Add(new Complex(0, -wc));
            }

            var tempPoles = new List<Complex>(poles.Count * 2);
            // First set of poles + conjugates
            tempPoles.AddRange(from pole in poles
                where !pole.Equals(Complex.Zero)
                let term1 = 0.5 * bw / pole
                let term2 = 0.5 * Complex.Sqrt((bw * bw) / (pole * pole) - (4 * wc * wc))
                select term1 + term2);
            // Second set of poles + conjugates
            tempPoles.AddRange(from pole in poles
                where !pole.Equals(Complex.Zero)
                let term1 = 0.5 * bw / pole
                let term2 = 0.5 * Complex.Sqrt((bw * bw) / (pole * pole) - (4 * wc * wc))
                select term1 - term2);

            // Copy converted zeros to output array
            zeros.Clear();
            foreach (var zero in ztmp)
                zeros.Add(zero);

            // Copy converted poles to output array
            poles.Clear();
            foreach (var pole in tempPoles)
                poles.Add(pole);

            return gain;
        }
    }
}