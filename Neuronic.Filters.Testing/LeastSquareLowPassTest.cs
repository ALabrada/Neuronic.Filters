using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Butterworth;
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
            const double cutoffFrequency = 25d;
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
            const double cutoffFrequency = 62d;
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
            const double cutoffFrequency = 500d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.006842283497823, 0.008914501335337, 0.014813041026146, 0.024103805930681, 0.035966867791804, 0.049278511611899,
                0.062729444628505, 0.074966038560548, 0.084738567584520, 0.091039394663403, 0.093215086738670, 0.091039394663403,
                0.084738567584520, 0.074966038560548, 0.062729444628505, 0.049278511611899, 0.035966867791804, 0.024103805930681,
                0.014813041026146, 0.008914501335337, 0.006842283497823
            };

            TestLowPass(order, fs, cutoffFrequency, expected, error, true);
        }
    }
}
