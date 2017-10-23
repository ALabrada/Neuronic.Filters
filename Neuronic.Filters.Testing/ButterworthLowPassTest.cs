using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthLowPassTest
    {
        [TestMethod]
        public void TestLowPass8()
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
    }
}
