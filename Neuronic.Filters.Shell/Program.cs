using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Neuronic.Filters.Butterwoth;

namespace Neuronic.Filters.Shell
{
    class Program
    {
        static void Main(string[] args)
        {
            const double noiseAmp = 1e-7d;
            Console.WriteLine("Noise amplitude: {0}.", noiseAmp);
            const double signalAmp = 0;
            Console.WriteLine("Signal amplitude: {0}.", signalAmp);
            const double signalFreq = 150;
            Console.WriteLine("Signal frequency: {0} Hz.", signalFreq);
            const double fs = 31250d;
            Console.WriteLine("Sampling frequency: {0} Hz.", fs);
            const double time = 0.15d;
            Console.WriteLine("Analysis time: {0} s.", time);
            const double lowCut = 5d;
            Console.WriteLine("Filter low cut: {0} Hz.", lowCut);
            const double highCut = 2000d;
            Console.WriteLine("Filter high cut: {0} Hz.", highCut);

            Console.WriteLine("Generating data...");
            var samples = new float[(int) Math.Ceiling(fs * time)];
            if (signalAmp > 0)
                GenerateSinusoid(signalAmp, signalFreq, fs, samples);
            if (noiseAmp > 0)
                GenerateNoise(noiseAmp, samples);

            Console.WriteLine("Designing filters...");
            var lowPass = new LowPassButtersworthCoefficients(5, fs, highCut).Calculate();
            Console.WriteLine("Low pass filter:");
            foreach (var biquad in lowPass)
                Console.WriteLine(biquad);
            var highPass = new HighPassButtersworthCoefficients(3, fs, highCut).Calculate();
            Console.WriteLine("High pass filter:");
            foreach (var biquad in highPass)
                Console.WriteLine(biquad);

            Console.WriteLine("Processing data...");
            var result = new float[samples.Length];
            Array.Copy(samples, result, samples.Length);
            lowPass.Filter(result, 0, result, 0, result.Length);
            highPass.Filter(result, 0, result, 0, result.Length);

            Console.WriteLine("Result:");
            var min = double.MaxValue;
            var max = double.MinValue;
            var sum = 0d;
            var sqsum = 0d;
            for (int i = 0; i < result.Length; i++)
            {
                var value = result[i];
                min = Math.Min(min, value);
                max = Math.Max(max, value);
                sum += value;
                sqsum += value * value;
            }
            var mean = sum / result.Length;
            var variance = (sqsum - result.Length * mean * mean) / (result.Length - 1);

            Console.WriteLine("Minimum: {0}.", min);
            Console.WriteLine("Maximum: {0}.", max);
            Console.WriteLine("Mean value: {0}.", mean);
            Console.WriteLine("Variance: {0}.", variance);
            Console.WriteLine("Standard deviation: {0}.", Math.Sqrt(variance));
        }

        static void GenerateSinusoid(double amplitude, double frequency, double fs, float[] samples)
        {
            var theta = 2 * Math.PI * frequency / fs;
            for (int i = 0; i < samples.Length; i++)
                samples[i] += (float) (amplitude * Math.Cos(theta * i));
        }

        static void GenerateNoise(double amplitude, float[] samples, Random rnd = null)
        {
            rnd = rnd ?? new Random();
            for (int i = 0; i < samples.Length; i++)
            {
                var u1 = rnd.NextDouble();
                var u2 = rnd.NextDouble();
                var stdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
                samples[i] += (float)(amplitude * stdNormal);
            }
        }
    }
}
