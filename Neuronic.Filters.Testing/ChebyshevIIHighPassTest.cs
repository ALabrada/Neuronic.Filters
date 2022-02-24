using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Chebyshev;
using System.Linq;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ChebyshevIIHighPassTest
    {
        [TestMethod]
        public void TestHighPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double error = 1e-4;

            var coeff = new HighPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.854998, -1.709830, 0.854998, 1.000000, -1.883222, 0.886929
1.000000, -1.998433, 1.000000, 1.000000, -1.898298, 0.903343
1.000000, -1.996492, 1.000000, 1.000000, -1.927393, 0.934394
1.000000, -1.995120, 1.000000, 1.000000, -1.967963, 0.976469").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 10d;
            const double error = 1e-4;

            var coeff = new HighPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.995933, -1.991866, 0.995933, 1.000000, -1.997892, 0.997893
1.000000, -1.999999, 1.000000, 1.000000, -1.998035, 0.998036
1.000000, -1.999999, 1.000000, 1.000000, -1.998311, 0.998313
1.000000, -1.999998, 1.000000, 1.000000, -1.998702, 0.998706
1.000000, -1.999997, 1.000000, 1.000000, -1.999182, 0.999186
1.000000, -1.999996, 1.000000, 1.000000, -1.999717, 0.999722").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new HighPassChebyshevIICoefficients(order, fs, cutoffFrequency, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.959945, -1.919886, 0.959945, 1.000000, -1.984108, 0.984176
1.000000, -1.999966, 1.000000, 1.000000, -1.984682, 0.984779
1.000000, -1.999910, 1.000000, 1.000000, -1.985811, 0.985964
1.000000, -1.999837, 1.000000, 1.000000, -1.987462, 0.987687
1.000000, -1.999758, 1.000000, 1.000000, -1.989580, 0.989884
1.000000, -1.999686, 1.000000, 1.000000, -1.992096, 0.992474
1.000000, -1.999630, 1.000000, 1.000000, -1.994925, 0.995359
1.000000, -1.999600, 1.000000, 1.000000, -1.997966, 0.998430").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}
