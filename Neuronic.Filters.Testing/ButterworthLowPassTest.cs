﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Accord.Diagnostics;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterworth;
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

            var coeff = new LowPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                new Biquad(2.16573939147328e-12,4.33147878294656e-12,2.16573939147328e-12,1,-1.96760498595716,0.972608207502736),
                new Biquad(1,2,1,1,-1.91904037532835,0.923920106890933),
                new Biquad(1,2,1,1,-1.88346019867072,0.888249457039896),
                new Biquad(1,2,1,1,-1.86474911008349,0.869490789938577),
            };

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

            var coeff = new LowPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.LowPass12).ToList();

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

            var coeff = new LowPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.LowPass16).ToList();
            var expectedGain = 3.444463173412853e-25;

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

            var coeff = new LowPassButterworthCoefficients(order, fs, cutoffFrequency);
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

        [TestMethod]
        public void TestLowPass2()
        {
            const int order = 2;
            const double fs = 44100d;
            const double cutoffFrequency = 2000d;
            const double error = 1e-4;

            var coeff = new LowPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                new Biquad(0.016819150107057118, 0.033638300214114236, 0.016819150107057118, 1, -1.6010923941836190, 0.66836899461184751),
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
