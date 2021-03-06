﻿using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class IIRFilteringTest
    {
        private static void TestFilter(int order, int fs, int cycles, IBiquadCoefficients coeff,
            IEnumerable<double> frequencies, IEnumerable<double> validFrequencies)
        {
            var samples = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var originalSignal = new Signal(samples);

            var coeffs = new List<Biquad>();
            var gain = coeff.Calculate(coeffs);
            var chain = new BiquadChain(coeffs, gain);
            chain.Filter(samples, 0, samples, 0, samples.Length);
            var filteredSignal = new Signal(samples, 0, samples.Length);

            Array.Clear(samples, 0, samples.Length);
            foreach (var frequency in validFrequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var expectedSignal = new Signal(samples);

            var filteredCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, 0);
            var originalCorrelation = Signal.CrossCorrelation(originalSignal, filteredSignal, 0);

            for (int i = 3; i <= order; i++)
            {
                var crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
                crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, -i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
            }

            Assert.AreEqual(1, filteredCorrelation, 0.1);
            Assert.AreNotEqual(1, originalCorrelation, 0.1);
            Assert.IsTrue(originalCorrelation < filteredCorrelation);
        }

        [TestMethod]
        public void TestLowPassButtersworthFiltering()
        {
            const int order = 16;
            const int fs = 44100;
            const double cutoffFrequency = 2000d;
            const int cycles = 10;
            double[] frequencies = { 770, 5830 };
            var validFrequencies = frequencies.TakeWhile(f => f < cutoffFrequency);

            var coeff = new LowPassButtersworthCoefficients(order, fs, cutoffFrequency);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestHighPassButtersworthFiltering()
        {
            const int order = 16;
            const int fs = 44100;
            const double cutoffFrequency = 1000d;
            const int cycles = 10;
            double[] frequencies = { 330, 1870 };
            var validFrequencies = frequencies.SkipWhile(f => f < cutoffFrequency);

            var coeff = new HighPassButtersworthCoefficients(order, fs, cutoffFrequency);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestBandPassButtersworthFiltering()
        {
            const int order = 32;
            const int fs = 44100;
            const double f1 = 1000d;
            const double f2 = 4000d;
            const int cycles = 10;
            double[] frequencies = { 330, 1870, 9790 };
            var validFrequencies = frequencies.SkipWhile(f => f < f1).TakeWhile(f => f < f2);

            var coeff = new BandPassButtersworthCoefficients(order, fs, f1, f2);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestBandStopButtersworthFiltering()
        {
            const int order = 32;
            const int fs = 44100;
            const double f1 = 500d;
            const double f2 = 8000d;
            const int cycles = 10;
            double[] frequencies = { 330, 770, 1870, 5830, 9790 };
            var validFrequencies = frequencies.TakeWhile(f => f < f1).Concat(frequencies.SkipWhile(f => f < f2));

            var coeff = new BandStopButtersworthCoefficients(order, fs, f1, f2);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }
    }
}
