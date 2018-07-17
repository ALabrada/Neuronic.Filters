using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class FIRFilteringTest
    {
        private static void TestFilter(int order, int fs, int cycles, FiniteImpulseResponseCoefficients coeff,
            IEnumerable<double> frequencies, IEnumerable<double> validFrequencies)
        {
            var samples = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var originalSignal = new Signal(samples);

            var chain = coeff.Calculate();
            chain.Filter(samples, 0, samples, 0, samples.Length, zeroPhase: true);
            var filteredSignal = new Signal(samples, 0, samples.Length - chain.PhaseShift);

            Array.Clear(samples, 0, samples.Length);
            foreach (var frequency in validFrequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var expectedSignal = new Signal(samples);

            var filteredCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, 0);
            var originalCorrelation = Signal.CrossCorrelation(originalSignal, filteredSignal, 0);

            for (int i = 1; i <= order; i++)
            {
                var crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
                crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, -i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
            }
            
            Assert.AreEqual(1, filteredCorrelation, 0.3);
            Assert.AreNotEqual(1, originalCorrelation, 0.3);
            Assert.IsTrue(originalCorrelation < filteredCorrelation);
        }

        [TestMethod]
        public void TestLowPassLeastSquareFiltering()
        {
            const int order = 65;
            const int fs = 44100;
            const double cutoffFrequency = 500d;
            const int cycles = 10;
            double[] frequencies = { 330, 770, 1870, 5830, 9790 };
            var validFrequencies = frequencies.TakeWhile(f => f < cutoffFrequency);

            var coeff = new LowPassLeastSquareCoefficients(order, fs, cutoffFrequency);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestHighPassLeastSquareFiltering()
        {
            const int order = 65;
            const int fs = 44100;
            const double cutoffFrequency = 8000d;
            const int cycles = 10;
            double[] frequencies = { 330, 770, 1870, 5830, 9790 };
            var validFrequencies = frequencies.SkipWhile(f => f < cutoffFrequency);

            var coeff = new HighPassLeastSquareCoefficients(order, fs, cutoffFrequency);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestBandPassLeastSquareFiltering()
        {
            const int order = 257;
            const int fs = 44100;
            const double f1 = 1000d;
            const double f2 = 4000d;
            const int cycles = 10;
            double[] frequencies = { 330, 770, 1870, 5830, 9790 };
            var validFrequencies = frequencies.SkipWhile(f => f < f1).TakeWhile(f => f < f2);

            var coeff = new BandPassLeastSquareCoefficients(order, fs, f1, f2);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }

        [TestMethod]
        public void TestBandStopLeastSquareFiltering()
        {
            const int order = 257;
            const int fs = 44100;
            const double f1 = 500d;
            const double f2 = 8000d;
            const int cycles = 10;
            double[] frequencies = { 330, 770, 1870, 5830, 9790 };
            var validFrequencies = frequencies.TakeWhile(f => f < f1).Concat(frequencies.SkipWhile(f => f < f2));

            var coeff = new BandStopLeastSquareCoefficients(order, fs, f1, f2);

            TestFilter(order, fs, cycles, coeff, frequencies, validFrequencies);
        }
    }
}
