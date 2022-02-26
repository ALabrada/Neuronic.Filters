using System;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.Chebyshev;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ChebyshevIIBandPassTest
    {
        [TestMethod]
        public void TestBandPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double br = 5;
            const double error = 1e-3;

            var coeff = new BandPassChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.003969, -0.007915, 0.003969, 1.000000, -1.993346, 0.998448
1.000000, -1.995422, 1.000000, 1.000000, -1.993423, 0.998458
1.000000, -1.994742, 1.000000, 1.000000, -1.993880, 0.999021
1.000000, -1.995107, 1.000000, 1.000000, -1.994035, 0.999035
1.000000, -1.994804, 1.000000, 1.000000, -1.994371, 0.999521
1.000000, -1.995048, 1.000000, 1.000000, -1.994534, 0.999528
1.000000, -1.994823, 1.000000, 1.000000, -1.994709, 0.999858
1.000000, -1.995030, 1.000000, 1.000000, -1.994865, 0.999861").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 1000d;
            const double br = 10;
            const double error = 1e-2;

            var coeff = new BandPassChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.003918, -0.007662, 0.003918, 1.000000, -1.954538, 0.993199
1.000000, -1.967016, 1.000000, 1.000000, -1.955330, 0.993263
1.000000, -1.959520, 1.000000, 1.000000, -1.956396, 0.995499
1.000000, -1.963525, 1.000000, 1.000000, -1.957994, 0.995586
1.000000, -1.960295, 1.000000, 1.000000, -1.958253, 0.997436
1.000000, -1.962813, 1.000000, 1.000000, -1.959900, 0.997488
1.000000, -1.960597, 1.000000, 1.000000, -1.959430, 0.998584
1.000000, -1.962528, 1.000000, 1.000000, -1.960953, 0.998611
1.000000, -1.960736, 1.000000, 1.000000, -1.960164, 0.999284
1.000000, -1.962395, 1.000000, 1.000000, -1.961580, 0.999297
1.000000, -1.960794, 1.000000, 1.000000, -1.960677, 0.999781
1.000000, -1.962339, 1.000000, 1.000000, -1.962034, 0.999784").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 60;
            const double br = 1;
            const double error = 1e-3;

            var coeff = new BandPassChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.003971, -0.007942, 0.003971, 1.000000, -1.998895, 0.999044
1.000000, -1.999896, 1.000000, 1.000000, -1.998921, 0.999064
1.000000, -1.999837, 1.000000, 1.000000, -1.999204, 0.999356
1.000000, -1.999870, 1.000000, 1.000000, -1.999244, 0.999383
1.000000, -1.999844, 1.000000, 1.000000, -1.999469, 0.999621
1.000000, -1.999864, 1.000000, 1.000000, -1.999498, 0.999637
1.000000, -1.999847, 1.000000, 1.000000, -1.999623, 0.999775
1.000000, -1.999862, 1.000000, 1.000000, -1.999644, 0.999784
1.000000, -1.999848, 1.000000, 1.000000, -1.999712, 0.999863
1.000000, -1.999861, 1.000000, 1.000000, -1.999728, 0.999868
1.000000, -1.999849, 1.000000, 1.000000, -1.999767, 0.999918
1.000000, -1.999860, 1.000000, 1.000000, -1.999780, 0.999921
1.000000, -1.999849, 1.000000, 1.000000, -1.999806, 0.999956
1.000000, -1.999859, 1.000000, 1.000000, -1.999816, 0.999957
1.000000, -1.999850, 1.000000, 1.000000, -1.999836, 0.999986
1.000000, -1.999859, 1.000000, 1.000000, -1.999846, 0.999987").ToList();

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandPassSinusoid()
        {
            const int order = 16;
            const int fs = 44100;
            const int targetFrequency = 3;
            const double br = 5;
            const int cycles = 10;
            double[] frequencies =
                {65.406, 130.81, 261.63, 523.25, 1046.5, 2093.0, 4186.0, 8372.0};

            var signal = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, signal);
            var im = new double[signal.Length];

            var coeff = new BandPassChebyshevIICoefficients(order, fs, frequencies[targetFrequency] - br, frequencies[targetFrequency] + br, 48);
            var chain = coeff.Calculate();
            chain.Filter(signal, 0, signal, 0, signal.Length);

            var count = signal.Length / 2;
            FourierTransform2.FFT(signal, im, FourierTransform.Direction.Forward);
            Helpers.CalculateEnergy(signal, im, count);

            var maxEnergy = signal.Take(count).Max();
            var step = fs / (2d * count);
            for (int i = 1; i < count - 1; i++)
            {
                var freq = i * step;
                if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= 0.05 * maxEnergy)
                {
                    var peak = frequencies.FirstOrDefault(x => Math.Abs(freq - x) <= 1);
                    Assert.AreEqual(frequencies[targetFrequency], peak);
                }
            }
        }
    }
}
