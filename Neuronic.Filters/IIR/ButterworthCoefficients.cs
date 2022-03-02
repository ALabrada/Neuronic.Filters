using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters.IIR
{
    /// <summary>
    /// Base class for the butterwoth filters.
    /// </summary>
    public abstract class ButterworthCoefficients : IBiquadCoefficients
    {
        private readonly List<Complex> _poles;
        private readonly List<Complex> _zeros;

        /// <summary>
        /// Creates a new instance of <see cref="ButtersworthCoefficients"/>.
        /// </summary>
        /// <param name="filterOrder">The order of the filter.</param>
        /// <param name="fs">The sampling frequency.</param>
        protected ButterworthCoefficients(int filterOrder, double fs)
        {
            FilterOrder = filterOrder;
            SamplingFrequency = fs;

            _zeros = new List<Complex>(2 * filterOrder);
            _poles = new List<Complex>(2 * filterOrder);
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }
        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        /// <summary>
        /// Gets the number of biquad filters.
        /// </summary>
        protected abstract int NumFilters { get; }

        /// <summary>
        /// Gets the poles.
        /// </summary>
        protected IList<Complex> Poles => _poles;

        /// <summary>
        /// Gets the zeroes.
        /// </summary>
        protected IList<Complex> Zeros => _zeros;

        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>The gain of the filter.</returns>
        public virtual double Calculate(IList<Biquad> coeffs)
        {
            // Init internal state based on filter design requirements
            _poles.Clear();
            _zeros.Clear();

            // Get zeros &poles of prototype analogue low pass.
            var tempPoles = PrototypeAnalogLowPass(FilterOrder);
            // Copy tmppole into poles
            _poles.AddRange(tempPoles);

            // Convert prototype to target filter type (LP/HP/BP/BS) - S-plane
            var gain = ConvertPoles();

            // SANITY CHECK: Ensure poles are in the left half of the S-plane
            for (var i = 0; i < Poles.Count; i++)
                Debug.Assert(Poles[i].Real <= 0, "Error: poles must be in the left half plane.");

            // Map zeros & poles from S-plane to Z-plane
            var ba = new double[2 * Math.Max(Poles.Count, Zeros.Count) + 5];
            var preBLTgain = gain;
            ConvertS2Z(Zeros, Poles, ref gain);

            //Split up Z-plane poles and zeros into SOS
            ConvertZp2SOS(Zeros, Poles, ref gain, ba);

            // correct the overall gain
            CorrectOverallGain(gain, preBLTgain, ba);

            // Init biquad chain with coefficients from SOS
            var overallGain = ba[0];
            var numFilters = NumFilters;

            coeffs.Clear();
            for (var i = 0; i < numFilters; i++)
                coeffs.Insert(0, new Biquad(1d, ba[4 * i + 1], ba[4 * i + 2], 1d, ba[4 * i + 3], ba[4 * i + 4]));

            return overallGain;
        }

        /// <summary>
        /// Creates a <see cref="BiquadChain"/> using the results from <see cref="Calculate(System.Collections.Generic.IList{Neuronic.Filters.Biquad})"/>.
        /// </summary>
        /// <returns>A biquad chain.</returns>
        public BiquadChain Calculate()
        {
            var coeffs = new List<Biquad>((FilterOrder + 1) / 2);
            var gain = Calculate(coeffs);
            return new TransposedDirectFormIIBiquadChain(coeffs, gain);
        }

        /// <summary>
        /// Converts the prototype low-pass filter to S-plane.
        /// </summary>
        /// <returns>The current gain.</returns>
        protected abstract double ConvertPoles();

        /// <summary>
        /// Performs final adjustments based on the filter type.
        /// </summary>
        protected virtual void CorrectOverallGain(double gain, double preBLTgain, double[] ba)
        {
        }

        /// <summary>
        ///     // Lowpass analogue prototype. Places Butterworth poles evenly around the S-plane unit circle.
        /// </summary>
        /// <param name="filterOrder">The filter order.</param>
        /// <returns></returns>
        /// <remarks>
        ///     Reference: MATLAB buttap(filterOrder)
        /// </remarks>
        private static List<Complex> PrototypeAnalogLowPass(int filterOrder)
        {
            var poles = new List<Complex>(filterOrder + 1);
            var count = (filterOrder + 1) / 2;
            for (var k = 0; k < count; k++)
            {
                var theta = (2 * k + 1) * Math.PI / (2 * filterOrder);
                var real = -Math.Sin(theta);
                var imag = Math.Cos(theta);
                poles.Add(new Complex(real, imag));
                poles.Add(new Complex(real, -imag)); // conjugate
            }
            return poles;
        }

        /// <summary>
        ///     Bilinears the transform.
        /// </summary>
        /// <returns>Z = (2 + S) / (2 - S) is the S-plane to Z-plane bilinear transform</returns>
        /// <remarks>
        ///     Reference: http://en.wikipedia.org/wiki/Bilinear_transform
        /// </remarks>
        private static double BilinearTransform(IList<Complex> values, int index)
        {
            Complex two = 2d;
            var s = values[index];
            values[index] = (two + s) / (two - s);
            return Complex.Abs(two - s);
        }

        /// <summary>
        ///     Convert poles & zeros from S-plane to Z-plane via Bilinear Tranform (BLT)
        /// </summary>
        /// <param name="zeros">The zeros.</param>
        /// <param name="numZeros">The number zeros.</param>
        /// <param name="poles">The poles.</param>
        /// <param name="numPoles">The number poles.</param>
        /// <param name="gain">The gain.</param>
        private static void ConvertS2Z(IList<Complex> zeros, IList<Complex> poles, ref double gain)
        {
            var numPoles = poles.Count;
            var numZeros = zeros.Count;
            // blt zeros
            for (var i = 0; i < numZeros; i++)
                gain /= BilinearTransform(zeros, i);
            // blt poles
            for (var i = 0; i < numPoles; i++)
                gain *= BilinearTransform(poles, i);
        }

        /// <summary>
        ///     Convert filter poles and zeros to second-order sections
        /// </summary>
        /// <param name="zeros">The zeros.</param>
        /// <param name="numZeros">The number zeros.</param>
        /// <param name="poles">The poles.</param>
        /// <param name="numPoles">The number poles.</param>
        /// <param name="gain">The gain.</param>
        /// <param name="ba">The ba.</param>
        /// <returns></returns>
        /// <remarks>
        ///     Reference: http://www.mathworks.com/help/signal/ref/zp2sos.html
        /// </remarks>
        private static int ConvertZp2SOS(IList<Complex> zeros, IList<Complex> poles, ref double gain,
            double[] ba)
        {
            var filterOrder = Math.Max(zeros.Count, poles.Count);
            // Copy
            var zerosTempVec = new List<Complex>(filterOrder);
            zerosTempVec.AddRange(zeros);
            // Add zeros at -1, so if S-plane degenerate case where
            // numZeros = 0 will map to -1 in Z-plane.
            if (zerosTempVec.Count < filterOrder)
                zerosTempVec.AddRange(Enumerable.Repeat(new Complex(-1, 0), filterOrder - zerosTempVec.Count));

            // Copy
            var polesTempVec = new List<Complex>(filterOrder);
            polesTempVec.AddRange(poles);
            if (polesTempVec.Count < filterOrder)
                polesTempVec.AddRange(Enumerable.Repeat(new Complex(0, 0), filterOrder - polesTempVec.Count));

            ba[0] = gain; // store gain

            var numSOS = 0;
            for (var i = 0; i < filterOrder - 1; i += 2, numSOS++)
            {
                ba[4 * numSOS + 1] = -(zerosTempVec[i] + zerosTempVec[i + 1]).Real;
                ba[4 * numSOS + 2] = (zerosTempVec[i] * zerosTempVec[i + 1]).Real;
                ba[4 * numSOS + 3] = -(polesTempVec[i] + polesTempVec[i + 1]).Real;
                ba[4 * numSOS + 4] = (polesTempVec[i] * polesTempVec[i + 1]).Real;
            }

            // Odd filter order thus one pair of poles/zeros remains
            if (filterOrder % 2 == 1)
            {
                ba[4 * numSOS + 1] = -zerosTempVec[filterOrder - 1].Real;
                ba[4 * numSOS + 2] = ba[4 * numSOS + 4] = 0;
                ba[4 * numSOS + 3] = -polesTempVec[filterOrder - 1].Real;
                numSOS++;
            }

            return 1 + 4 * numSOS;
        }

        /// <summary>
        ///     Convert analog lowpass prototype poles to lowpass
        /// </summary>
        protected static double ConvertToLowPass(double freq, IList<Complex> poles, IList<Complex> zeros)
        {
            var gain = Math.Pow(freq, poles.Count);

            zeros.Clear();
            for (var i = 0; i < poles.Count; i++)
                poles[i] = freq * poles[i];

            return gain;
        }

        /// <summary>
        /// Convert lowpass poles & zeros to highpass with Wc = f2, use:  hp_S = Wc / lp_S;
        /// </summary>
        /// <param name="freq">The freq.</param>
        /// <param name="poles">The poles.</param>
        /// <param name="zeros">The zeros.</param>
        /// <returns></returns>
        protected static double ConvertToHighPass(double freq, IList<Complex> poles, IList<Complex> zeros)
        {
            // Calculate gain
            var prodz = zeros.Aggregate(new Complex(1, 0), (current, zero) => current * -zero);
            var prodp = poles.Aggregate(new Complex(1, 0), (current, pole) => current * -pole);
            var gain = prodz.Real / prodp.Real;

            // Convert LP poles to HP
            for (int i = 0; i < poles.Count; i++)
                if (!poles[i].Equals(Complex.Zero))
                    poles[i] = freq / poles[i];

            // Init with zeros, no non-zero values to convert
            zeros.Clear();
            for (int i = 0; i < poles.Count; i++)
                zeros.Add(Complex.Zero);

            return gain;
        }

        protected static double ConvertBandPass(double f1, double f2, IList<Complex> poles,
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

        protected static double ConvertBandStop(double f1, double f2, IList<Complex> poles,
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