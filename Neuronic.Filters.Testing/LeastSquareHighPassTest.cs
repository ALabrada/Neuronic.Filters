using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class LeastSquareHighPassTest
    {
        private static void TestHighPass(int order, double fs, double cutoffFrequency, double[] expected, double error, bool scale)
        {
            var coeff = new HighPassLeastSquareCoefficients(order, fs, cutoffFrequency) {UseScaling = scale};
            var chain = coeff.Calculate();

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000090689927722, -0.000190284900112, -0.000451030080377, -0.000773333542649, -0.001034085222857, 0.998775626423610,
                -0.001034085222857, -0.000773333542649, -0.000451030080377, -0.000190284900112, -0.000090689927722
            };

            TestHighPass(order, fs, cutoffFrequency, expected, error, true);
        }

        [TestMethod]
        public void TestHighPass15()
        {
            const int order = 15;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000090698899030, -0.000130401280119, -0.000243462595537, -0.000412675045610, -0.000612279711861, -0.000811888541729,
                -0.000981111426445, -0.001094183748515, 0.998956808033798, -0.001094183748515, -0.000981111426445, -0.000811888541729,
                -0.000612279711861, -0.000412675045610, -0.000243462595537, -0.000130401280119, -0.000090698899030
            };

            TestHighPass(order, fs, cutoffFrequency, expected, error, true);
        }

        [TestMethod]
        public void TestHighPass20NoScale()
        {
            const int order = 20;
            const double fs = 44100d;
            const double cutoffFrequency = 1000d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000887971477626, -0.001142485081560, -0.001877436030273, -0.003025327634321, -0.004476536684647, -0.006090135201398,
                -0.007707962549592, -0.009170520265426, -0.010333067030265, -0.011080264393001, 0.988662131519274, -0.011080264393001,
                -0.010333067030265, -0.009170520265426, -0.007707962549592, -0.006090135201398, -0.004476536684647, -0.003025327634321,
                -0.001877436030273, -0.001142485081560, -0.000887971477626
            };

            TestHighPass(order, fs, cutoffFrequency, expected, error, false);
        }

        [TestMethod]
        public void TestHighPass25NoScale()
        {
            const int order = 25;
            const double fs = 44100d;
            const double cutoffFrequency = 1000d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000874962181838, -0.001026640486894, -0.001466227032755, -0.002172034352309, -0.003105911455549, -0.004215443429900,
                -0.005437054289052, -0.006699824424834, -0.007929787394688, -0.009054439000098, -0.010007176814770, -0.010731391539375,
                -0.011183952651993, 0.988662131519274, -0.011183952651993, -0.010731391539375, -0.010007176814770, -0.009054439000098,
                -0.007929787394688, -0.006699824424834, -0.005437054289052, -0.004215443429900, -0.003105911455549, -0.002172034352309,
                -0.001466227032755, -0.001026640486894, -0.000874962181838 
            };

            TestHighPass(order, fs, cutoffFrequency, expected, error, false);
        }
    }
}