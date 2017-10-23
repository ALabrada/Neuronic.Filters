using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthHighPassTest
    {
        [TestMethod]
        public void TestHighPass8()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double error = 1e-4;

            var coeff = new HighPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.HighPass08).Reverse().ToList();
            var expectedGain = 0.833143245502442;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 10d;
            const double error = 1e-4;

            var coeff = new HighPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.HighPass12).Reverse().ToList();
            var expectedGain = 0.992506754917111;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new HighPassButtersworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.HighPass16).Reverse().ToList();
            var expectedGain = 0.902520827102739;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}