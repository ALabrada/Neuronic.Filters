using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters
{
    /// <summary>
    /// Uses a sequence of biquad filters to implement a higher-order filter.
    /// </summary>
    public class BiquadChain : IZeroPhaseDigitalFilter
#if NET40
        , IEnumerable<Biquad>
#else
        , IReadOnlyList<Biquad>
#endif
    {
        private readonly Biquad[] _coeffs;
        private readonly double[] _yn;
        private readonly double[] _yn1;
        private readonly double[] _yn2;
        private double _xn1;
        private double _xn2;

        /// <summary>
        /// Initializes a new instance of the <see cref="BiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        public BiquadChain(IList<Biquad> coefficients, double gain)
        {
            _coeffs = new Biquad[coefficients.Count];
            coefficients.CopyTo(_coeffs, 0);
            Gain = gain;
            _yn = new double[_coeffs.Length];
            _yn1 = new double[_coeffs.Length];
            _yn2 = new double[_coeffs.Length];
        }

        /// <summary>
        /// Gets the overall gain of the filter.
        /// </summary>
        public double Gain { get; }

        /// <summary>
        /// Reset's the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        public void Reset()
        {
            _xn1 = 0;
            _xn2 = 0;
            Array.Clear(_yn, 0, Count);
            Array.Clear(_yn1, 0, Count);
            Array.Clear(_yn2, 0, Count);
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(float* input, float* output, int count, int stride)
        {
            fixed (Biquad* coeffs = _coeffs)
            {
                fixed (double* yn = _yn, yn1 = _yn1, yn2 = _yn2)
                {
                    for (int n = 0; n < count; n++)
                    {
                        var xn = *input * Gain;

                        yn[0] = coeffs[0].B0 * xn + coeffs[0].B1 * _xn1 + coeffs[0].B2 * _xn2
                                + coeffs[0].A1 * yn1[0] + coeffs[0].A2 * yn2[0];

                        for (int i = 1; i < Count; i++)
                            yn[i] = coeffs[i].B0 * yn[i - 1] + coeffs[i].B1 * yn1[i - 1] + coeffs[i].B2 * yn2[i - 1]
                                    + coeffs[i].A1 * yn1[i] + coeffs[i].A2 * yn2[i];

                        // Shift delay line elements.
                        for (int i = 0; i < Count; i++)
                        {
                            yn2[i] = yn1[i];
                            yn1[i] = yn[i];
                        }
                        _xn2 = _xn1;
                        _xn1 = xn;

                        // Store result and stride
                        *output = (float)yn[Count - 1];

                        input += stride;
                        output += stride;
                    }
                }
            }
        }

        /// <summary>
        /// Executes two sweeps of the filter (forward and backward). This is a zero-phase filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <remarks>
        /// This method will reset the filter.
        /// </remarks>
        public unsafe void Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count, int stride=1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
                Reset();
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
                Reset();
            }
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(float[] input, int inputIndex, float[] output, int outputIndex, int count,
            int stride=1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
            }
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(double* input, double* output, int count, int stride)
        {
            fixed (Biquad* coeffs = _coeffs)
            {
                fixed (double* yn = _yn, yn1 = _yn1, yn2 = _yn2)
                {
                    for (int n = 0; n < count; n++)
                    {
                        var xn = *input * Gain;

                        yn[0] = coeffs[0].B0 * xn + coeffs[0].B1 * _xn1 + coeffs[0].B2 * _xn2
                                + coeffs[0].A1 * yn1[0] + coeffs[0].A2 * yn2[0];

                        for (int i = 1; i < Count; i++)
                            yn[i] = coeffs[i].B0 * yn[i - 1] + coeffs[i].B1 * yn1[i - 1] + coeffs[i].B2 * yn2[i - 1]
                                    + coeffs[i].A1 * yn1[i] + coeffs[i].A2 * yn2[i];

                        // Shift delay line elements.
                        for (int i = 0; i < Count; i++)
                        {
                            yn2[i] = yn1[i];
                            yn1[i] = yn[i];
                        }
                        _xn2 = _xn1;
                        _xn1 = xn;

                        // Store result and stride
                        *output = yn[Count - 1];

                        input += stride;
                        output += stride;
                    }
                }
            }
        }

        /// <summary>
        /// Executes two sweeps of the filter (forward and backward). This is a zero-phase filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <remarks>
        /// This method will reset the filter.
        /// </remarks>
        public unsafe void Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count, int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
                Reset();
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
                Reset();
            }
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(double[] input, int inputIndex, double[] output, int outputIndex, int count,
            int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
            }
        }

        /// <summary>Returns an enumerator that iterates through the collection.</summary>
        /// <returns>An enumerator that can be used to iterate through the collection.</returns>
        /// <filterpriority>1</filterpriority>
        public IEnumerator<Biquad> GetEnumerator()
        {
            return _coeffs.Cast<Biquad>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        /// <summary>Gets the number of elements in the collection.</summary>
        /// <returns>The number of elements in the collection. </returns>
        public int Count => _coeffs.Length;

        /// <summary>Gets the element at the specified index in the read-only list.</summary>
        /// <returns>The element at the specified index in the read-only list.</returns>
        /// <param name="index">The zero-based index of the element to get. </param>
        public Biquad this[int index] => _coeffs[index];
    }
}