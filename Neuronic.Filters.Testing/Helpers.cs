using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Neuronic.Filters.Testing
{
    public static class Helpers
    {
        internal static void ValidateBiquad(Biquad expected, Biquad actual, double error)
        {
            Assert.AreEqual(expected.B0, actual.B0, error);
            Assert.AreEqual(expected.B1, actual.B1, error);
            Assert.AreEqual(expected.B2, actual.B2, error);
            Assert.AreEqual(expected.A1, actual.A1, error);
            Assert.AreEqual(expected.A2, actual.A2, error);
        }

        internal static IEnumerable<Biquad> LoadScript(string script)
        {
            foreach (var line in script.Split(new char[] { '\n', '\r'}, StringSplitOptions.RemoveEmptyEntries))
            {
                var values = new List<double>(6);
                values.AddRange(line.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => double.Parse(s, NumberStyles.Any, CultureInfo.CurrentCulture)));
                yield return new Biquad(values[0], values[1], values[2], values[3], values[4], values[5]);
            }
        }

        internal static IEnumerable<Biquad> LoadCsv(string csv)
        {
            foreach (var line in csv.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries))
            {
                var values = new List<double>(6);
                values.AddRange(line.Split(new char[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => double.Parse(s, NumberStyles.Any, CultureInfo.InvariantCulture)));
                yield return new Biquad(values[0], values[1], values[2], values[3], values[4], values[5]);
            }
        }

        internal static void GenerateSinusoid(double frequency, double fs, double[] samples)
        {
            var theta = 2 * Math.PI * frequency / fs;
            for (int i = 0; i < samples.Length; i++)
                samples[i] += Math.Cos(theta * i);
        }

        internal static void CalculateEnergy(double[] re, double[] im, int count)
        {
            for (int i = 0; i < count; i++)
                re[i] = 2 * Math.Sqrt(re[i] * re[i] + im[i] * im[i]);
        }

        internal static IEnumerable<double> GetTransitionBands(this BiquadChain chain, double fs, double db = -3.0)
        {
            var len = (int)Math.Floor(fs / 2.0);
            var response = new double[len];
            for (int i = 0; i < response.Length; i++)
            {
                var impulse = chain.GetResponse(i / fs);
                response[i] = 20 * Math.Log10(impulse.Magnitude + 1e-16);
            }

            for (int i = 1; i < response.Length; i++)
            {
                var prev = response[i - 1] - db;
                var current = response[i] - db;
                if (prev < 0 && current >= 0 ||
                    prev > 0 && current <= 0)
                    yield return i;
            }
        }
    }
}