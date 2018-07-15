using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class WindowBasedHighPassTest
    {
        [TestMethod]
        public void TestHighPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new HighPassWindowBasedCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                -0.000090689927722, -0.000190284900112, -0.000451030080377, -0.000773333542649, -0.001034085222857, 0.998775626423610,
                -0.001034085222857, -0.000773333542649, -0.000451030080377, -0.000190284900112, -0.000090689927722
            };

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }
    }
}