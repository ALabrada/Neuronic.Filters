using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Accord.Diagnostics;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthLowPassTest
    {
        [TestMethod]
        public void TestLowPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double error = 1e-4;

            var coeff = new LowPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                new Biquad(1, 2, 1, 1, -1.96762058043629, 0.97261960500367),
                new Biquad(1, 2, 1, 1, -1.91907529383978, 0.92395098208778),
                new Biquad(1, 2, 1, 1, -1.88350864178159, 0.88829396780773),
                new Biquad(1, 2, 1, 1, -1.86480445083537, 0.86954225616013),
            };
            var expectedGain = 2.158589092625337e-12;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 2000d;
            const double error = 1e-2;

            var coeff = new LowPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.LowPass12).Reverse().ToList();
            var expectedGain = 7.343036284781901e-10;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 300d;
            const double error = 1e-4;

            var coeff = new LowPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.LowPass16).Reverse().ToList();
            var expectedGain = 3.444463173412853e-25;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPassSinusoid()
        {
            const int order = 16;
            const int fs = 44100;
            const double cutoffFrequency = 400d;
            const int cycles = 10;
            double[] frequencies =
                {65.406, 130.81, 261.63, 523.25, 1046.5, 2093.0, 4186.0, 8372.0};

            var signal = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, signal);
            var im = new double[signal.Length];

            var coeff = new LowPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();
            chain.Filter(signal, 0, signal, 0, signal.Length);

            var count = signal.Length / 2;
            FourierTransform2.FFT(signal, im, FourierTransform.Direction.Forward);
            Helpers.CalculateEnergy(signal, im, count);

            var maxEnergy = signal.Take(count).Max();
            var step = fs / (2d * count);
            var peakSet = new HashSet<double>();
            for (int i = 1; i < count - 1; i++)
            {
                var freq = i * step;
                if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= 0.001 * maxEnergy)
                {
                    var peak = frequencies.FirstOrDefault(x => Math.Abs(freq - x) <= 1);
                    Assert.AreNotEqual(0, peak);
                    peakSet.Add(peak);
                }
            }
            Assert.IsTrue(peakSet.SetEquals(frequencies.Where(x => x < cutoffFrequency)));
        }
    }
}
