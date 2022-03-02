using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.IIR;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ChebyshevIIBandStopTest
    {
        [TestMethod]
        public void TestBandStop08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double br = 5;
            const double error = 1e-3;

            var coeff = new BandStopChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.996871, -1.988667, 0.996871, 1.000000, -1.993703, 0.998799
1.000000, -1.994948, 1.000000, 1.000000, -1.993761, 0.998805
1.000000, -1.994871, 1.000000, 1.000000, -1.993832, 0.998977
1.000000, -1.994984, 1.000000, 1.000000, -1.993995, 0.998991
1.000000, -1.994843, 1.000000, 1.000000, -1.994131, 0.999314
1.000000, -1.995011, 1.000000, 1.000000, -1.994368, 0.999328
1.000000, -1.994827, 1.000000, 1.000000, -1.994554, 0.999758
1.000000, -1.995026, 1.000000, 1.000000, -1.994822, 0.999765").ToList();
            var expectedGain = 0.996359732766342;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStop12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 1000d;
            const double br = 10;
            const double error = 1e-2;

            var coeff = new BandStopChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.991883, -1.945553, 0.991883, 1.000000, -1.959391, 0.997890
1.000000, -1.961674, 1.000000, 1.000000, -1.959624, 0.997896
1.000000, -1.961280, 1.000000, 1.000000, -1.959306, 0.998028
1.000000, -1.961866, 1.000000, 1.000000, -1.959988, 0.998045
1.000000, -1.961105, 1.000000, 1.000000, -1.959375, 0.998302
1.000000, -1.962038, 1.000000, 1.000000, -1.960458, 0.998325
1.000000, -1.960962, 1.000000, 1.000000, -1.959597, 0.998694
1.000000, -1.962177, 1.000000, 1.000000, -1.961001, 0.998717
1.000000, -1.960860, 1.000000, 1.000000, -1.959955, 0.999178
1.000000, -1.962276, 1.000000, 1.000000, -1.961581, 0.999195
1.000000, -1.960807, 1.000000, 1.000000, -1.960426, 0.999719
1.000000, -1.962326, 1.000000, 1.000000, -1.962159, 0.999725").ToList();
            var expectedGain = 0.985211119221826;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStop16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 60;
            const double br = 1;
            const double error = 1e-3;

            var coeff = new BandStopChebyshevIICoefficients(order, fs, cutoffFrequency - br, cutoffFrequency + br, 48);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(
                @"0.999183, -1.998220, 0.999183, 1.000000, -1.999694, 0.999840
1.000000, -1.999855, 1.000000, 1.000000, -1.999696, 0.999841
1.000000, -1.999853, 1.000000, 1.000000, -1.999699, 0.999846
1.000000, -1.999856, 1.000000, 1.000000, -1.999703, 0.999847
1.000000, -1.999852, 1.000000, 1.000000, -1.999710, 0.999857
1.000000, -1.999857, 1.000000, 1.000000, -1.999717, 0.999860
1.000000, -1.999851, 1.000000, 1.000000, -1.999726, 0.999875
1.000000, -1.999858, 1.000000, 1.000000, -1.999735, 0.999878
1.000000, -1.999851, 1.000000, 1.000000, -1.999747, 0.999897
1.000000, -1.999858, 1.000000, 1.000000, -1.999758, 0.999900
1.000000, -1.999850, 1.000000, 1.000000, -1.999773, 0.999923
1.000000, -1.999859, 1.000000, 1.000000, -1.999785, 0.999926
1.000000, -1.999850, 1.000000, 1.000000, -1.999802, 0.999953
1.000000, -1.999859, 1.000000, 1.000000, -1.999814, 0.999954
1.000000, -1.999850, 1.000000, 1.000000, -1.999833, 0.999984
1.000000, -1.999859, 1.000000, 1.000000, -1.999844, 0.999985").ToList();
            var expectedGain = 0.997950883359409;

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestBandStopSinusoid()
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

            var coeff = new BandStopChebyshevIICoefficients(order, fs, frequencies[targetFrequency] - br, frequencies[targetFrequency] + br, 48);
            var chain = coeff.Calculate();
            chain.Filter(signal, 0, signal, 0, signal.Length);

            var count = signal.Length / 2;
            FourierTransform2.FFT(signal, im, FourierTransform.Direction.Forward);
            Helpers.CalculateEnergy(signal, im, count);

            var maxEnergy = signal.Take(count).Max();
            var step = fs / (2d * count);
            var peakSet = new HashSet<double>();
            for (int i = 1; i < count - 1; i++)
            {
                var freq = i * step;
                if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= 0.01 * maxEnergy)
                {
                    var peak = frequencies.FirstOrDefault(x => Math.Abs(freq - x) <= 1);
                    Assert.AreNotEqual(0, peak);
                    peakSet.Add(peak);
                }
            }
            Assert.IsTrue(peakSet.SetEquals(frequencies.Except(Enumerable.Repeat(frequencies[targetFrequency], 1))));
        }
    }
}
