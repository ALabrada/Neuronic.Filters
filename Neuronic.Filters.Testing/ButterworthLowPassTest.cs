using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Accord.Diagnostics;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.IIR;
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
                new Biquad(0.000000, 0.000000, 0.000000, 1.000000, -1.864749, 0.869491),
                new Biquad(1.000000, 2.000000, 1.000000, 1.000000, -1.883460, 0.888249),
                new Biquad(1.000000, 2.000000, 1.000000, 1.000000, -1.919040, 0.923920),
                new Biquad(1.000000, 2.000000, 1.000000, 1.000000, -1.967605, 0.972608),
            };

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 2.0);

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

            var expected = Helpers.LoadCsv(Resources.LowPass12).ToList();

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 1.0);

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

            var expected = Helpers.LoadCsv(Resources.LowPass16).ToList();
            var expectedGain = 3.444463173412853e-25;

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 1.0);

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
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

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 5.0);

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
