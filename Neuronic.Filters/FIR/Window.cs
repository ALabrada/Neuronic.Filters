using System;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters.FIR
{
    public interface IWindow
    {
        void CopyTo(double[] array, int start, int count);
        void ApplyTo(IList<double> coeffs);
    }

    public class DiscreteWindow : IWindow
    {
        private readonly double[] _samples;

        public DiscreteWindow(IList<double> samples)
        {
            _samples = samples.ToArray();
        }

        public double this[int index] => _samples[index];

        public void CopyTo(double[] array, int start, int count)
        {
            if (count != _samples.Length)
                throw new ArgumentException(nameof(count));
            _samples.CopyTo(array, start);
        }

        public void ApplyTo(IList<double> coeffs)
        {
            if (coeffs.Count != _samples.Length)
                throw new ArgumentException(nameof(coeffs));
            var b = coeffs;
            for (int i = 0; i < b.Count; i++)
                b[i] *= this[i];
        }
    }

    public abstract class ContinuousWindow : IWindow
    {
        public abstract double this[double t] { get; }

        public void CopyTo(double[] array, int start, int count)
        {
            if (count < 2)
                return;
            var resolution = 1d / (count - 1);
            for (int i = 0; i < count; i++)
                array[start + i] = this[resolution * i];
        }

        public void ApplyTo(IList<double> coeffs)
        {
            var resolution = 1d / (coeffs.Count - 1);
            var b = coeffs;
            for (int i = 0; i < b.Count; i++)
                b[i] *= this[resolution * i];
        }
    }

    public class Window : ContinuousWindow
    {
        private readonly Func<double, double> _func;

        public Window(Func<double, double> func)
        {
            if (func == null) throw new ArgumentNullException(nameof(func));
            _func = func;
        }

        public readonly static Window Rect = new Window(t => 1d);

        public readonly static Window Hamming = new Window(t => 0.54 - 0.46 * Math.Cos(2 * Math.PI * t));

        public readonly static Window Blackman =
            new Window(t => 0.42 - 0.5 * Math.Cos(Math.PI * t) + 0.08 * Math.Cos(2 * Math.PI * t));

        public readonly static Window FlatTop =
            new Window(t => 0.21557895 -
                            0.41663158 * Math.Cos(2 * Math.PI * t) +
                            0.277263158 * Math.Cos(4 * Math.PI * t) -
                            0.083578947 * Math.Cos(6 * Math.PI * t) +
                            0.006947368 + Math.Cos(8 * Math.PI * t));

        public readonly static Window Hann = new Window(t => 0.5 * (1 - Math.Cos(2 * Math.PI * t)));

        public override double this[double t]
        {
            get { return _func(t); }
        }
    }
}