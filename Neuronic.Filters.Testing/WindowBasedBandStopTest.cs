using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class WindowBasedBandStopTest
    {
        [TestMethod]
        public void TestBandStop10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double f1 = 200d;
            const double f2 = 800d;
            const double error = 1e-4;

            var coeff = new BandStopWindowBasedCoefficients(order, fs, f1, f2);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                -0.000562766728964, -0.001182663056506, -0.002806702624320, -0.004816583398478, -0.006444026035493, 1.031625483687525,
                -0.006444026035493, -0.004816583398478, -0.002806702624320, -0.001182663056506, -0.000562766728964
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }
    }
}