using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class WindowBasedHighPassTest
    {
        [TestMethod]
        public void TestHighPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new HighPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                -0.000090689927722, -0.000190284900112, -0.000451030080377, -0.000773333542649, -0.001034085222857, 0.998775626423610,
                -0.001034085222857, -0.000773333542649, -0.000451030080377, -0.000190284900112, -0.000090689927722
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPassSinusoid()
        {
            const int order = 32;
            const int fs = 44100;
            const double cutoffFrequency = 4000d;
            const int cycles = 10;
            double[] frequencies = { 330, 1870, 5830, 9790 };

            var samples = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var originalSignal = new Signal(samples);

            var coeff = new HighPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();
            chain.Filter(samples, 0, samples, 0, samples.Length, zeroPhase: true);
            var filteredSignal = new Signal(samples, 0, samples.Length - chain.PhaseShift);

            Array.Clear(samples, 0, samples.Length);
            foreach (var frequency in frequencies.SkipWhile(f => f < cutoffFrequency))
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var expectedSignal = new Signal(samples);

            var filteredCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, 0);
            Assert.AreEqual(1, filteredCorrelation, 0.3);

            var originalCorrelation = Signal.CrossCorrelation(originalSignal, filteredSignal, 0);
            Assert.AreNotEqual(1, originalCorrelation, 0.3);
            Assert.IsTrue(originalCorrelation < filteredCorrelation);

            for (int i = 1; i <= order; i++)
            {
                var crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
                crossCorrelation = Signal.CrossCorrelation(expectedSignal, filteredSignal, -i);
                Assert.IsTrue(crossCorrelation < filteredCorrelation);
            }
        }
    }
}