﻿using System;
using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    public abstract class FourierSeriesCoefficients
    {
        public FourierSeriesCoefficients(int n, double fs)
        {
            if (n <= 0)
                throw new ArgumentOutOfRangeException(nameof(n));
            if (fs <= 0)
                throw new ArgumentOutOfRangeException(nameof(fs));

            FilterOrder = n;
            SamplingFrequency = fs;
        }

        /// <summary>
        /// Gets the order of the filter.
        /// </summary>
        public int FilterOrder { get; }
        /// <summary>
        /// Gets the sampling frequency of the filter.
        /// </summary>
        public double SamplingFrequency { get; set; }

        public TapChain Calculate()
        {
            var coeffs = new List<double>(FilterOrder);
            Calculate(coeffs);
            return new TapChain(coeffs);
        }

        public virtual void Calculate(IList<double> coeffs)
        {
            coeffs.Clear();
            for (var n = 0; n <= FilterOrder; n++)
                coeffs.Add(CalculateTap(n));
        }

        protected abstract double CalculateTap(int n);
    }

    public class LowPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda;

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

        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? _lambda / Math.PI
                : Math.Sin(mm * _lambda) / (mm * Math.PI);
        }
    }

    public class HighPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda;

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

        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? 1d - _lambda / Math.PI
                : -Math.Sin(mm * _lambda) / (mm * Math.PI);
        }
    }

    public class BandPassFourierSeriesCoefficients : FourierSeriesCoefficients
    {
        private readonly double _lambda, _phi;

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

        protected override double CalculateTap(int n)
        {
            var mm = n - FilterOrder / 2d;
            return mm == 0
                ? (_phi - _lambda) / Math.PI
                : (Math.Sin(mm * _phi) - Math.Sin(mm * _lambda)) / (mm * Math.PI);
        }
    }
}