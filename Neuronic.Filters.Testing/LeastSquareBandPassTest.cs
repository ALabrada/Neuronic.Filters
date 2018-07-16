using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class LeastSquareBandPassTest
    {
        [TestMethod]
        public void TestBandPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double f1 = 200d;
            const double f2 = 800d;
            const double error = 1e-4;

            var coeff = new BandPassLeastSquareCoefficients(order, fs, f1, f2);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                0.014555018853135, 0.030587599085740, 0.072590662364310, 0.124572862190256, 0.166663940153847, 0.182748392063313,
                0.166663940153847, 0.124572862190256, 0.072590662364310, 0.030587599085740, 0.014555018853135
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }
    }
}