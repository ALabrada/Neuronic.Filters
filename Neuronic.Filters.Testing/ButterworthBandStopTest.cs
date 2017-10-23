using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthBandStopTest
    {
        [TestMethod]
        public void TestBandStop8()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double br = 5;
            const double error = 1e-2;

            var coeff = new BandStopButtersworthCoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadScript(Resources.BandStop08).Reverse().ToList();
            var expectedGain = 0.996359732766342;

            Assert.AreEqual(expectedGain, chain.Gain, error);
            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
