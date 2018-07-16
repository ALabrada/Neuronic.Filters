using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public partial class WindowBasedLowPassTest
    {
        [TestMethod]
        public void TestLowPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new LowPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                0.014597902574489, 0.030629205502635, 0.072600048724931, 0.124479619696451, 0.166451509181935,
                0.18248342863912, 0.166451509181935, 0.124479619696451, 0.072600048724931, 0.030629205502635,
                0.014597902574489
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass15()
        {
            const int order = 15;
            const double fs = 44100d;
            const double cutoffFrequency = 250d;
            const double error = 1e-4;

            var coeff = new LowPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                0.009773912665721, 0.014635366052850, 0.028378523543657, 0.048630291918410, 0.071890208924709,
                0.094136103024869, 0.111520234113506, 0.121035359756277, 0.121035359756277, 0.111520234113506,
                0.094136103024869, 0.071890208924709, 0.048630291918410, 0.028378523543657, 0.014635366052850,
                0.009773912665721
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass20()
        {
            const int order = 20;
            const double fs = 44100d;
            const double cutoffFrequency = 1250d;
            const double error = 1e-4;

            var coeff = new LowPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                0.007151722151675, 0.009222810505808, 0.015186933196523, 0.024516688062819, 0.036333673358802, 0.049495542655145,
                0.062711289192240, 0.074672841958675, 0.084189241153656, 0.090309309752530, 0.092419896024252, 0.090309309752530,
                0.084189241153656, 0.074672841958675, 0.062711289192240, 0.049495542655145, 0.036333673358802, 0.024516688062819,
                0.015186933196523, 0.009222810505808, 0.007151722151675
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
            double[] frequencies = {330, 1870, 5830, 9790 };

            var samples = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, samples);
            var originalSignal = new Signal(samples);

            var coeff = new LowPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();
            chain.Filter(samples, 0, samples, 0, samples.Length, zeroPhase: true);
            var filteredSignal = new Signal(samples, 0, samples.Length - chain.PhaseShift);

            Array.Clear(samples, 0, samples.Length);
            foreach (var frequency in frequencies.TakeWhile(f => f < cutoffFrequency))
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
