using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Chebyshev;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ChebyshevIILowPassTest
    {

        [TestMethod]
        public void TestLowPass2()
        {
            const int order = 2;
            const double fs = 44100d;
            const double cutoffFrequency = 2000d;
            const double error = 1e-4;

            var coeff = new LowPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                new Biquad(1, -1.8418888082730389, 1, 1, -1.9638716093794057, 0.96451523313991561),
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
