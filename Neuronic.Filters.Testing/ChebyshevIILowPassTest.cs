using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.IIR;
using System.Linq;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ChebyshevIILowPassTest
    {
        [TestMethod]
        public void TestLowPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double error = 1e-4;

            var coeff = new LowPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.003542, -0.006627, 0.003542, 1.000000, -1.850372, 0.856450
1.000000, -1.983612, 1.000000, 1.000000, -1.902676, 0.907306
1.000000, -1.992667, 1.000000, 1.000000, -1.950079, 0.953551
1.000000, -1.994727, 1.000000, 1.000000, -1.983090, 0.986058
").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 2000d;
            const double error = 1e-2;

            var coeff = new LowPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.004224, 0.003363, 0.004224, 1.000000, -0.961310, 0.242000
1.000000, -1.149174, 1.000000, 1.000000, -1.179991, 0.413183
1.000000, -1.614136, 1.000000, 1.000000, -1.428074, 0.609010
1.000000, -1.763423, 1.000000, 1.000000, -1.615363, 0.760065
1.000000, -1.822796, 1.000000, 1.000000, -1.745316, 0.870017
1.000000, -1.845222, 1.000000, 1.000000, -1.840403, 0.958141").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestLowPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 300d;
            const double error = 1e-4;

            var coeff = new LowPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.003088, -0.005107, 0.003088, 1.000000, -1.733557, 0.752389
1.000000, -1.957258, 1.000000, 1.000000, -1.814016, 0.827649
1.000000, -1.983684, 1.000000, 1.000000, -1.885766, 0.894795
1.000000, -1.990975, 1.000000, 1.000000, -1.929719, 0.935987
1.000000, -1.993917, 1.000000, 1.000000, -1.955793, 0.960505
1.000000, -1.995325, 1.000000, 1.000000, -1.972269, 0.976106
1.000000, -1.996028, 1.000000, 1.000000, -1.983739, 0.987102
1.000000, -1.996328, 1.000000, 1.000000, -1.992751, 0.995910").ToList();
            var expectedGain = 3.444463173412853e-25;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
