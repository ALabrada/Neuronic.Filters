using System;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters.FIR
{
    /// <summary>
    /// Abstraction of a windowing function.
    /// </summary>
    public interface IWindow
    {
        /// <summary>
        /// Copies the window's samples to an array.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="start">The start index.</param>
        /// <param name="count">The count.</param>
        void CopyTo(double[] array, int start, int count);
        /// <summary>
        /// Applies the windowing function to a tap vector.
        /// </summary>
        /// <param name="coeffs">The FIR coefficients (taps).</param>
        void ApplyTo(IList<double> coeffs);
    }

    /// <summary>
    /// A discrete windowing function.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.IWindow" />
    public class DiscreteWindow : IWindow
    {
        private readonly double[] _samples;

        /// <summary>
        /// Initializes a new instance of the <see cref="DiscreteWindow"/> class.
        /// </summary>
        /// <param name="samples">The samples.</param>
        public DiscreteWindow(IList<double> samples)
        {
            _samples = samples.ToArray();
        }

        /// <summary>
        /// Gets the sample at the specified index.
        /// </summary>
        /// <value>
        /// The sample.
        /// </value>
        /// <param name="index">The index.</param>
        /// <returns>The window's sample at <paramref name="index"/>.</returns>
        public double this[int index] => _samples[index];

        /// <summary>
        /// Copies the window's samples to an array.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="start">The start index.</param>
        /// <param name="count">The count.</param>
        /// <exception cref="ArgumentException">Thrown when <paramref name="count"/>
        /// does not match the size of the window function.</exception>
        public void CopyTo(double[] array, int start, int count)
        {
            if (count != _samples.Length)
                throw new ArgumentException(nameof(count));
            _samples.CopyTo(array, start);
        }

        /// <summary>
        /// Applies the windowing function to a tap vector.
        /// </summary>
        /// <param name="coeffs">The FIR coefficients (taps).</param>
        /// <exception cref="ArgumentException">Throw when the tap vector's size 
        /// does not match the size of the window function.</exception>
        public void ApplyTo(IList<double> coeffs)
        {
            if (coeffs.Count != _samples.Length)
                throw new ArgumentException(nameof(coeffs));
            var b = coeffs;
            for (int i = 0; i < b.Count; i++)
                b[i] *= this[i];
        }
    }

    /// <summary>
    /// A continuous windowing function defined in [0,1].
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.IWindow" />
    public abstract class ContinuousWindow : IWindow
    {
        /// <summary>
        /// Gets or sets the function's value at the specified <paramref name="t"/> time instant.
        /// </summary>
        /// <param name="t">The time instant.</param>
        /// <returns>Function's value at <paramref name="t"/>.</returns>
        public abstract double this[double t] { get; }

        /// <summary>
        /// Copies the window's samples to an array.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="start">The start index.</param>
        /// <param name="count">The count.</param>
        public void CopyTo(double[] array, int start, int count)
        {
            if (count < 2)
                return;
            var resolution = 1d / (count - 1);
            for (int i = 0; i < count; i++)
                array[start + i] = this[resolution * i];
        }

        /// <summary>
        /// Applies the windowing function to a tap vector.
        /// </summary>
        /// <param name="coeffs">The FIR coefficients (taps).</param>
        public void ApplyTo(IList<double> coeffs)
        {
            var resolution = 1d / (coeffs.Count - 1);
            var b = coeffs;
            for (int i = 0; i < b.Count; i++)
                b[i] *= this[resolution * i];
        }
    }

    /// <summary>
    /// A convenient way of defining continuous windows using delegates.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.ContinuousWindow" />
    public class Window : ContinuousWindow
    {
        private readonly Func<double, double> _func;

        /// <summary>
        /// Initializes a new instance of the <see cref="Window"/> class.
        /// </summary>
        /// <param name="func">The function.</param>
        /// <exception cref="ArgumentNullException">func</exception>
        public Window(Func<double, double> func)
        {
            if (func == null) throw new ArgumentNullException(nameof(func));
            _func = func;
        }

        /// <summary>
        /// A rectangular window.
        /// </summary>
        public readonly static Window Rect = new Window(t => 1d);

        /// <summary>
        /// A Hamming window.
        /// </summary>
        public readonly static Window Hamming = new Window(t => 0.54 - 0.46 * Math.Cos(2 * Math.PI * t));

        /// <summary>
        /// A Backman window.
        /// </summary>
        public readonly static Window Blackman =
            new Window(t => 0.42 - 0.5 * Math.Cos(Math.PI * t) + 0.08 * Math.Cos(2 * Math.PI * t));

        /// <summary>
        /// A Flattop window.
        /// </summary>
        public readonly static Window FlatTop =
            new Window(t => 0.21557895 -
                            0.41663158 * Math.Cos(2 * Math.PI * t) +
                            0.277263158 * Math.Cos(4 * Math.PI * t) -
                            0.083578947 * Math.Cos(6 * Math.PI * t) +
                            0.006947368 + Math.Cos(8 * Math.PI * t));

        /// <summary>
        /// A Hann window.
        /// </summary>
        public readonly static Window Hann = new Window(t => 0.5 * (1 - Math.Cos(2 * Math.PI * t)));

        /// <summary>
        /// Gets the <see cref="System.Double"/> with the specified t.
        /// </summary>
        /// <value>
        /// The <see cref="System.Double"/>.
        /// </value>
        /// <param name="t">The t.</param>
        /// <returns></returns>
        public override double this[double t] => _func(t);
    }
}