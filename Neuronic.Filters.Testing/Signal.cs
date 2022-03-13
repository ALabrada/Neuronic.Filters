using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Neuronic.Filters.Testing
{
    class Signal : IReadOnlyList<double>
    {
        private readonly double[] _samples;

        public Signal(IList<double> samples, int start, int count)
        {
            _samples = new double[count];

            var sum = 0d;
            var sqSum = 0d;
            for (int i = 0; i < count; i++)
            {
                var x = _samples[i] = samples[start + i];
                sum += x;
                sqSum += x * x;
            }

            Mean = sum / count;
            Variance = (sqSum - count * Mean * Mean) / (count - 1);
            StandardDeviation = Math.Sqrt(Variance);
        }

        public Signal(IList<double> samples) : this(samples, 0, samples.Count)
        {
        }

        public double Mean { get; }

        public double Variance { get; }

        public double StandardDeviation { get; }

        public int Count => _samples.Length;

        public double this[int index] => _samples[index];

        public IEnumerator<double> GetEnumerator()
        {
            return _samples.Cast<double>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public static double CrossCorrelation(Signal x, Signal y, int m)
        {
            var sum = 0d;
            var count = 0;
            for (int i = 0; i < x.Count; i++)
                if (i + m >= 0 && i + m < y.Count)
                {
                    sum += (x[i] - x.Mean) * (y[i + m] - y.Mean);
                    count++;
                }
            var productMean = sum / count;
            return productMean / (x.StandardDeviation * y.StandardDeviation);
        }

        public string Export()
        {
            return string.Join(",", this.Select(x => x.ToString(CultureInfo.InvariantCulture)));
        }

        public Signal Extract(int start, int count)
        {
            if (start < 0 || start >= Count)
                throw new IndexOutOfRangeException(nameof(start));
            if (count < 0 || start + count > Count)
                throw new ArgumentOutOfRangeException(nameof(count));
            return new Signal(_samples, start, count);
        }

        public override string ToString()
        {
            return $"Avg: {Mean}, Var: {Variance}, Std: {StandardDeviation}";
        }
    }
}
