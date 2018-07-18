using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class LeastSquareBandStopTest
    {
        private static void TestBandStop(int order, double fs, double f1, double f2, double[] expected, double error, bool scale)
        {
            var coeff = new BandStopLeastSquareCoefficients(order, fs, f1, f2) {UseScaling = scale};
            var chain = coeff.Calculate();

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStop10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double f1 = 50d;
            const double f2 = 200d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000562766728964, -0.001182663056506, -0.002806702624320, -0.004816583398478, -0.006444026035493, 1.031625483687525,
                -0.006444026035493, -0.004816583398478, -0.002806702624320, -0.001182663056506, -0.000562766728964
            };

            TestBandStop(order, fs, f1, f2, expected, error, true);
        }

        [TestMethod]
        public void TestBandStop15()
        {
            const int order = 15;
            const double fs = 44100d;
            const double f1 = 50d;
            const double f2 = 200d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000571907913556, -0.000824433303171, -0.001542768881151, -0.002620099192127, -0.003893558824616, -0.005169249245692,
                -0.006252169415976, -0.006976400131888, 1.055701173816355, -0.006976400131888, -0.006252169415976, -0.005169249245692,
                -0.003893558824616, -0.002620099192127, -0.001542768881151, -0.000824433303171, -0.000571907913556
            };

            TestBandStop(order, fs, f1, f2, expected, error, true);
        }

        [TestMethod]
        public void TestBandStop20NoScale()
        {
            const int order = 20;
            const double fs = 44100d;
            const double f1 = 500d;
            const double f2 = 2000d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.000931236795395, 0.000190622994486, -0.001466502528377, -0.005311811837723, -0.012154001120775, -0.022023919244603,
                -0.034059793940273, -0.046636512751245, -0.057714392914079, -0.065319410870160, 0.931972789115646, -0.065319410870160,
                -0.057714392914079, -0.046636512751245, -0.034059793940273, -0.022023919244603, -0.012154001120775, -0.005311811837723,
                -0.001466502528377, 0.000190622994486, 0.000931236795395
            };

            TestBandStop(order, fs, f1, f2, expected, error, false);
        }

        [TestMethod]
        public void TestBandStop25NoScale()
        {
            const int order = 25;
            const double fs = 44100d;
            const double f1 = 500d;
            const double f2 = 2000d;
            const double error = 1e-6;

            var expected = new[]
            {
                0.002610762679355, 0.002547817016009, 0.002682879569589, 0.002277864053854, 0.000518219582751, -0.003292766490307,
                -0.009546274924833, -0.018190328663607, -0.028676702804130, -0.040009577706685, -0.050891314300402, -0.059939197723603,
                -0.065930664875598, 0.931972789115646, -0.065930664875598, -0.059939197723603, -0.050891314300402, -0.040009577706685,
                -0.028676702804130, -0.018190328663607, -0.009546274924833, -0.003292766490307, 0.000518219582751, 0.002277864053854,
                0.002682879569589, 0.002547817016009, 0.002610762679355
            };

            TestBandStop(order, fs, f1, f2, expected, error, false);
        }
    }
}