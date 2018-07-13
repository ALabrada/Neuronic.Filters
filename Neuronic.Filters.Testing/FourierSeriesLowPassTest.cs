using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.FIR;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class FourierSeriesLowPassTest
    {
        [TestMethod]
        public void TestLowPassSinusoid()
        {
            const int order = 32;
            const int fs = 44100;
            const double cutoffFrequency = 400d;
            const int cycles = 10;
            double[] frequencies =
                {65.406, 130.81, 261.63, 523.25, 1046.5, 2093.0, 4186.0, 8372.0};

            var signal = new double[cycles * fs];
            foreach (var frequency in frequencies)
                Helpers.GenerateSinusoid(frequency, fs, signal);
            var im = new double[signal.Length];

            var coeff = new LowPassFourierSeriesCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();
            var filteredSignal = new double[signal.Length - chain.PhaseShift];
            chain.Filter(signal, 0, filteredSignal, 0, signal.Length, zeroPhase: true);

            var count = filteredSignal.Length / 2;
            FourierTransform2.FFT(filteredSignal, im, FourierTransform.Direction.Forward);
            Helpers.CalculateEnergy(filteredSignal, im, count);

            var maxEnergy = filteredSignal.Take(count).Max();
            var step = fs / (2d * count);
            var peakSet = new HashSet<double>();
            for (int i = 1; i < count - 1; i++)
            {
                var freq = i * step;
                if (filteredSignal[i] > filteredSignal[i - 1] && filteredSignal[i] > filteredSignal[i + 1] &&
                    filteredSignal[i] >= 0.001 * maxEnergy)
                {
                    var peak = frequencies.FirstOrDefault(x => Math.Abs(freq - x) <= 1);
                    Assert.AreNotEqual(0, peak);
                    peakSet.Add(peak);
                }
            }
            Assert.IsTrue(peakSet.SetEquals(frequencies.Where(x => x < cutoffFrequency)));
        }
    }
}
