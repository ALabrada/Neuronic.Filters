﻿using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class FourierSeriesBandPassTest
    {
        private static void TestBandPass(int order, double fs, double f1, double f2, double[] expected, double error, bool scale, IWindow window)
        {
            var coeff = new BandPassFourierSeriesCoefficients(order, fs, f1, f2) { Window = window, UseScaling = scale };
            var chain = coeff.Calculate();

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass50Hamming()
        {
            const int order = 50;
            const double fs = 200d;
            const double f1 = 10d;
            const double f2 = 30d;
            const double error = 1e-3;

            var expected = new[]
            {
                -0.002030, -0.001701, -0.000651, 0.000588, 0.001034, -0.000000, -0.001709, -0.001580, 0.002738, 0.010499, 0.016827, 0.015822, 0.006236, -0.005463, -0.009029, 0.000000, 0.012968, 0.011328, -0.018878, -0.071211, -0.115734, -0.115087, -0.051159, 0.056783, 0.158026, 0.199304, 0.158026, 0.056783, -0.051159, -0.115087, -0.115734, -0.071211, -0.018878, 0.011328, 0.012968, 0.000000, -0.009029, -0.005463, 0.006236, 0.015822, 0.016827, 0.010499, 0.002738, -0.001580, -0.001709, -0.000000, 0.001034, 0.000588, -0.000651, -0.001701, -0.002030
            };

            TestBandPass(order, fs, f1, f2, expected, error, true, Window.Hamming);
        }

        [TestMethod]
        public void TestBandPass50Rect()
        {
            const int order = 50;
            const double fs = 200d;
            const double f1 = 10d;
            const double f2 = 30d;
            const double error = 1e-5;

            var expected = new[]
            {
                -0.025089, -0.020108, -0.006818, 0.005178, 0.007467, -0.000000, -0.008253, -0.006329, 0.009224, 0.030162, 0.041815, 0.034471, 0.012062, -0.009494, -0.014255, 0.000000, 0.017423, 0.014241, -0.022401, -0.080433, -0.125444, -0.120649, -0.052268, 0.056963, 0.156805, 0.197047, 0.156805, 0.056963, -0.052268, -0.120649, -0.125444, -0.080433, -0.022401, 0.014241, 0.017423, 0.000000, -0.014255, -0.009494, 0.012062, 0.034471, 0.041815, 0.030162, 0.009224, -0.006329, -0.008253, -0.000000, 0.007467, 0.005178, -0.006818, -0.020108, -0.025089
            };

            TestBandPass(order, fs, f1, f2, expected, error, true, Window.Rect);
        }

        [TestMethod]
        public void TestBandPass99HammingNoScale()
        {
            const int order = 99;
            const double fs = 200d;
            const double f1 = 10d;
            const double f2 = 30d;
            const double error = 1e-3;

            var expected = new[]
            {
                0.000153, 0.000283, -0.000000, -0.000633, -0.001245, -0.001384, -0.000866, 0.000000, 0.000565, 0.000358, -0.000407, -0.000828, 0.000000, 0.002082, 0.004210, 0.004735, 0.002961, -0.000000, -0.001882, -0.001167, 0.001295, 0.002572, -0.000000, -0.006153, -0.012150, -0.013359, -0.008179, 0.000000, 0.005011, 0.003061, -0.003354, -0.006591, 0.000000, 0.015559, 0.030664, 0.033760, 0.020776, -0.000000, -0.013041, -0.008133, 0.009161, 0.018680, -0.000000, -0.049315, -0.105712, -0.130414, -0.094184, -0.000000, 0.113018, 0.189386, 0.189386, 0.113018, -0.000000, -0.094184, -0.130414, -0.105712, -0.049315, -0.000000, 0.018680, 0.009161, -0.008133, -0.013041, -0.000000, 0.020776, 0.033760, 0.030664, 0.015559, 0.000000, -0.006591, -0.003354, 0.003061, 0.005011, 0.000000, -0.008179, -0.013359, -0.012150, -0.006153, -0.000000, 0.002572, 0.001295, -0.001167, -0.001882, -0.000000, 0.002961, 0.004735, 0.004210, 0.002082, 0.000000, -0.000828, -0.000407, 0.000358, 0.000565, 0.000000, -0.000866, -0.001384, -0.001245, -0.000633, -0.000000, 0.000283, 0.000153
            };

            TestBandPass(order, fs, f1, f2, expected, error, false, Window.Hamming);
        }

        [TestMethod]
        public void TestBandPass99RectNoScale()
        {
            const int order = 99;
            const double fs = 200d;
            const double f1 = 10d;
            const double f2 = 30d;
            const double error = 1e-5;

            var expected = new[]
            {
                0.001913, 0.003503, -0.000000, -0.007170, -0.013143, -0.013438, -0.007665, 0.000000, 0.004094, 0.002339, -0.002398, -0.004413, 0.000000, 0.009135, 0.016845, 0.017334, 0.009953, -0.000000, -0.005393, -0.003105, 0.003211, 0.005961, -0.000000, -0.012582, -0.023451, -0.024408, -0.014188, 0.000000, 0.007901, 0.004620, -0.004857, -0.009183, 0.000000, 0.020207, 0.038581, 0.041242, 0.024697, -0.000000, -0.014772, -0.009020, 0.009970, 0.019986, -0.000000, -0.051294, -0.108729, -0.132890, -0.095260, -0.000000, 0.113254, 0.189430, 0.189430, 0.113254, -0.000000, -0.095260, -0.132890, -0.108729, -0.051294, -0.000000, 0.019986, 0.009970, -0.009020, -0.014772, -0.000000, 0.024697, 0.041242, 0.038581, 0.020207, 0.000000, -0.009183, -0.004857, 0.004620, 0.007901, 0.000000, -0.014188, -0.024408, -0.023451, -0.012582, -0.000000, 0.005961, 0.003211, -0.003105, -0.005393, -0.000000, 0.009953, 0.017334, 0.016845, 0.009135, 0.000000, -0.004413, -0.002398, 0.002339, 0.004094, 0.000000, -0.007665, -0.013438, -0.013143, -0.007170, -0.000000, 0.003503, 0.001913
            };

            TestBandPass(order, fs, f1, f2, expected, error, false, Window.Rect);
        }
    }
}