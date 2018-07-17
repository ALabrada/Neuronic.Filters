using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class LeastSquareBandPassTest
    {
        private static void TestBandPass(int order, double fs, double f1, double f2, double[] expected, double error, bool scale)
        {
            var coeff = new BandPassLeastSquareCoefficients(order, fs, f1, f2) {UseScaling = scale};
            var chain = coeff.Calculate();

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass10()
        {
            const int order = 10;
            const double fs = 44100d;
            const double f1 = 200d;
            const double f2 = 800d;
            const double error = 1e-6;
            
            var expected = new[]
            {
                0.014555018853135, 0.030587599085740, 0.072590662364310, 0.124572862190256, 0.166663940153847, 0.182748392063313,
                0.166663940153847, 0.124572862190256, 0.072590662364310, 0.030587599085740, 0.014555018853135
            };

            TestBandPass(order, fs, f1, f2, expected, error, true);
        }

        [TestMethod]
        public void TestHighPass15()
        {
            const int order = 15;
            const double fs = 44100d;
            const double f1 = 200d;
            const double f2 = 800d;
            const double error = 1e-6;
            
            var expected = new[]
            {
                0.009713480273588, 0.014578561315120, 0.028324409676720, 0.048617631390348, 0.071966240485030, 0.094328752308441,
                0.111822003642067, 0.121402793552637, 0.121402793552637, 0.111822003642067, 0.094328752308441, 0.071966240485030,
                0.048617631390348, 0.028324409676720, 0.014578561315120, 0.009713480273588
            };

            TestBandPass(order, fs, f1, f2, expected, error, true);
        }

        [TestMethod]
        public void TestHighPass20NoScale()
        {
            const int order = 20;
            const double fs = 44100d;
            const double f1 = 2000d;
            const double f2 = 8000d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.000931236795395, -0.000190622994486, 0.001466502528377, 0.005311811837723, 0.012154001120775, 0.022023919244603,
                0.034059793940273, 0.046636512751245, 0.057714392914079, 0.065319410870160, 0.068027210884354, 0.065319410870160,
                0.057714392914079, 0.046636512751245, 0.034059793940273, 0.022023919244603, 0.012154001120775, 0.005311811837723,
                0.001466502528377, -0.000190622994486, -0.000931236795395
            };

            TestBandPass(order, fs, f1, f2, expected, error, false);
        }

        [TestMethod]
        public void TestHighPass25NoScale()
        {
            const int order = 25;
            const double fs = 44100d;
            const double f1 = 2000d;
            const double f2 = 8000d;
            const double error = 1e-6;

            var expected = new[]
            {
                -0.002414893105356, -0.002262878856118, -0.002204256302043, -0.001407698158938, 0.000987275252323, 0.005650591595040,
                0.012865155759272, 0.022401613010615, 0.033498140963964, 0.044956292636187, 0.055337727921033, 0.063223179040726,
                0.067479791657420, 0.067479791657420, 0.063223179040726, 0.055337727921033, 0.044956292636187, 0.033498140963964,
                0.022401613010615, 0.012865155759272, 0.005650591595040, 0.000987275252323, -0.001407698158938, -0.002204256302043,
                -0.002262878856118, -0.002414893105356
            };

            TestBandPass(order, fs, f1, f2, expected, error, false);
        }
    }
}