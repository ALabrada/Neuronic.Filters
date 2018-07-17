﻿using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterwoth;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class LeastSquareLowPassTest
    {
        private static void TestLowPass(int order, double fs, double cutoffFrequency, double[] expected, double error, bool scale)
        {
            var coeff = new LowPassLeastSquareCoefficients(order, fs, cutoffFrequency) {UseScaling = scale};
            var chain = coeff.Calculate();

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.014597902574489, 0.030629205502635, 0.072600048724931, 0.124479619696451, 0.166451509181935,
                0.18248342863912, 0.166451509181935, 0.124479619696451, 0.072600048724931, 0.030629205502635,
                0.014597902574489
            };

            TestLowPass(order, fs, cutoffFrequency, expected, error, true);
        }

        [TestMethod]
        public void TestLowPass15()
        {
            const int order = 15;
            const double fs = 44100d;
            const double cutoffFrequency = 250d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.009773912665721, 0.014635366052850, 0.028378523543657, 0.048630291918410, 0.071890208924709,
                0.094136103024869, 0.111520234113506, 0.121035359756277, 0.121035359756277, 0.111520234113506,
                0.094136103024869, 0.071890208924709, 0.048630291918410, 0.028378523543657, 0.014635366052850,
                0.009773912665721
            };

            TestLowPass(order, fs, cutoffFrequency, expected, error, true);
        }

        [TestMethod]
        public void TestLowPass20()
        {
            const int order = 20;
            const double fs = 44100d;
            const double cutoffFrequency = 1250d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.007151722151675, 0.009222810505808, 0.015186933196523, 0.024516688062819, 0.036333673358802, 0.049495542655145,
                0.062711289192240, 0.074672841958675, 0.084189241153656, 0.090309309752530, 0.092419896024252, 0.090309309752530,
                0.084189241153656, 0.074672841958675, 0.062711289192240, 0.049495542655145, 0.036333673358802, 0.024516688062819,
                0.015186933196523, 0.009222810505808, 0.007151722151675
            };

            TestLowPass(order, fs, cutoffFrequency, expected, error, true);
        }
    }
}
