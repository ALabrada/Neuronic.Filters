using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterworth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthBandStopTest
    {
        [TestMethod]
        public void TestBandStop08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double br = 5;
            const double error = 1e-3;

            var coeff = new BandStopButterworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandStop08).Reverse().ToList();
            var expectedGain = 0.996359732766342;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStop12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 1000d;
            const double br = 10;
            const double error = 1e-2;

            var coeff = new BandStopButterworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandStop12).Reverse().ToList();
            var expectedGain = 0.985211119221826;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStop16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 60;
            const double br = 1;
            const double error = 1e-3;

            var coeff = new BandStopButterworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandStop16).Reverse().ToList();
            var expectedGain = 0.997950883359409;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStopSinusoid()
        {
            const int order = 16;
            const int fs = 44100;
            const int targetFrequency = 3;
            const double br = 5;
            const int cycles = 10;
            double[] frequencies =
                {65.406, 130.81, 261.63, 523.25, 1046.5, 2093.0, 4186.0, 8372.0};

            var signal = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, signal);
            var im = new double[signal.Length];

            var coeff = new BandStopButterworthCoefficients(order, fs, frequencies[targetFrequency] - br, frequencies[targetFrequency] + br);
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
                if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= 0.01 * maxEnergy)
                {
                    var peak = frequencies.FirstOrDefault(x => Math.Abs(freq - x) <= 1);
                    Assert.AreNotEqual(0, peak);
                    peakSet.Add(peak);
                }
            }
            Assert.IsTrue(peakSet.SetEquals(frequencies.Except(Enumerable.Repeat(frequencies[targetFrequency], 1))));
        }
    }
}
