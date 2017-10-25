using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthBandPassTest
    {
        [TestMethod]
        public void TestBandPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double br = 5;
            const double error = 1e-3;

            var coeff = new BandPassButtersworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandPass08).Reverse().ToList();
            var expectedGain = 6.541867072214145e-26;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 1000d;
            const double br = 10;
            const double error = 1e-2;

            var coeff = new BandPassButtersworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandPass12).Reverse().ToList();
            var expectedGain = 2.883388964711794e-33;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 60;
            const double br = 1;
            const double error = 1e-3;

            var coeff = new BandPassButtersworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandPass16).Reverse().ToList();
            var expectedGain = 7.114329431453330e-60;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
