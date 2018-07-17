using System;
using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    /// <summary>
    /// Abstraction of a FIR filter designer that uses the Fourier Series method.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.WindowBasedCoefficients" />
    public abstract class FourierSeriesCoefficients : WindowBasedCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="FourierSeriesCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        public FourierSeriesCoefficients(int n, double fs) : base(n, fs)
        {
        }

        /// <summary>
        /// Calculates coefficients (taps) for a FIR filter.
        /// </summary>
        /// <param name="coeffs">The collection that will hold the coefficients.</param>
        public override void Calculate(IList<double> coeffs)
        {
            coeffs.Clear();
            for (var n = 0; n <= FilterOrder; n++)
                coeffs.Add(CalculateTap(n));

            ApplyWindow(coeffs);
        }

        /// <summary>
        /// Calculates the tap at the specified index.
        /// </summary>
        /// <param name="n">The index.</param>
        protected abstract double CalculateTap(int n);
    }

    /// <summary>
    /// Low-pass FIR filter designer that uses the Fourier Series method.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.FourierSeriesCoefficients" />
    public class LowPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda;

        /// <summary>
        /// Initializes a new instance of the <see cref="LowPassFourierSeriesCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="fx">The cutoff frequency.</param>
        /// <exception cref="ArgumentOutOfRangeException">fx</exception>
        public LowPassFourierSeriesCoefficients(int n, double fs, double fx) : base(n, fs)
        {
            if (fx <= 0 || fx >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(fx));

            CutoffFrequency = fx;
            _lambda = Math.PI * fx / (fs / 2);
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <summary>
        /// Calculates the tap at the specified index.
        /// </summary>
        /// <param name="n">The index.</param>
        /// <returns></returns>
        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? _lambda / Math.PI
                : Math.Sin(mm * _lambda) / (mm * Math.PI);
        }
    }

    /// <summary>
    /// High-pass FIR filter designer that uses the Fourier Series method.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.FourierSeriesCoefficients" />
    public class HighPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda;

        /// <summary>
        /// Initializes a new instance of the <see cref="HighPassFourierSeriesCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="fx">The cutoff frequency.</param>
        /// <exception cref="ArgumentOutOfRangeException">fx</exception>
        public HighPassFourierSeriesCoefficients(int n, double fs, double fx) : base(n, fs)
        {
            if (fx <= 0 || fx >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(fx));

            CutoffFrequency = fx;
            _lambda = Math.PI * fx / (fs / 2);
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }

        /// <summary>
        /// Calculates the tap at the specified index.
        /// </summary>
        /// <param name="n">The index.</param>
        /// <returns></returns>
        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? 1d - _lambda / Math.PI
                : -Math.Sin(mm * _lambda) / (mm * Math.PI);
        }
    }

    /// <summary>
    /// Band-pass FIR filter designer that uses the Fourier Series method.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.FourierSeriesCoefficients" />
    public class BandPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda, _phi;

        /// <summary>
        /// Initializes a new instance of the <see cref="BandPassFourierSeriesCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The first lower frequency.</param>
        /// <param name="f2">The first upper frequency.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// f1
        /// or
        /// f2
        /// </exception>
        public BandPassFourierSeriesCoefficients(int n, double fs, double f1, double f2) : base(n, fs)
        {
            if (f1 <= 0 || f1 >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(f1));
            if (f2 <= 0 || f2>= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(f2));

            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;

            _lambda = Math.PI * f1 / (fs / 2);
            _phi = Math.PI * f2 / (fs / 2);
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }

        /// <summary>
        /// Calculates the tap at the specified index.
        /// </summary>
        /// <param name="n">The index.</param>
        /// <returns></returns>
        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? (_phi - _lambda) / Math.PI
                : (Math.Sin(mm * _phi) - Math.Sin(mm * _lambda)) / (mm * Math.PI);
        }
    }

    /// <summary>
    /// Band-stop FIR filter designer that uses the Fourier Series method.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.FourierSeriesCoefficients" />
    public class BandStopFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda, _phi;

        /// <summary>
        /// Initializes a new instance of the <see cref="BandStopFourierSeriesCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        /// <param name="f1">The first lower frequency.</param>
        /// <param name="f2">The first upper frequency.</param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// f1
        /// or
        /// f2
        /// </exception>
        public BandStopFourierSeriesCoefficients(int n, double fs, double f1, double f2) : base(n, fs)
        {
            if (f1 <= 0 || f1 >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(f1));
            if (f2 <= 0 || f2 >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(f2));

            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;

            _lambda = Math.PI * f1 / (fs / 2);
            _phi = Math.PI * f2 / (fs / 2);
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }

        /// <summary>
        /// Calculates the tap at the specified index.
        /// </summary>
        /// <param name="n">The index.</param>
        /// <returns></returns>
        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? 1d - (_phi - _lambda) / Math.PI
                : (Math.Sin(mm * _lambda) - Math.Sin(mm * _phi)) / (mm * Math.PI);
        }
    }
}