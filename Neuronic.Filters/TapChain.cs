using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters
{
    /// <summary>
    /// Represents a linear-phase finite impulse response (FIR) digital filter.
    /// </summary>
    public class TapChain : IZeroPhaseDigitalFilter
#if NET40
        , IEnumerable<double>
#else
        , IReadOnlyList<double>
#endif
    {
        private readonly double[] _taps;
        private readonly double[] _sr;

        /// <summary>
        /// Initializes a new instance of the <see cref="TapChain"/> class.
        /// </summary>
        /// <param name="taps">The filter coefficients or taps.</param>
        public TapChain(IList<double> taps)
        {
            _taps = taps.ToArray();
            _sr = new double[_taps.Length];
        }

        /// <summary>
        /// Gets the number of elements in the collection.
        /// </summary>
        public int Count => _taps.Length;

        /// <summary>
        /// Gets the group delay of the linear-phase filter.
        /// </summary>
        /// <value>
        /// The group delay of the linear-phase filter.
        /// </value>
        public int PhaseShift => Count / 2;

        /// <summary>
        /// Gets the <see cref="System.Double"/> at the specified index.
        /// </summary>
        /// <value>
        /// The <see cref="System.Double"/>.
        /// </value>
        /// <param name="index">The index.</param>
        /// <returns></returns>
        public double this[int index]
        {
            get { return _taps[index]; }
        }

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>
        /// An enumerator that can be used to iterate through the collection.
        /// </returns>
        public IEnumerator<double> GetEnumerator()
        {
            return _taps.Cast<double>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        /// <summary>
        /// Reset's the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        public void Reset()
        {
            Array.Clear(_sr, 0, Count);
        }

        private unsafe double DoSample(double* taps, double* sr, double dataSample)
        {
            var tapCount = Count;

            // Shift right
            for (int i = tapCount - 1; i >= 1; i--)
                sr[i] = sr[i - 1];
            sr[0] = dataSample;

            // Convolution
            var result = 0d;
            for (var i = 0; i < tapCount; i++)
                result += sr[i] * taps[i];
            return result;
        }

        /// <summary>
        /// Filters the specified single-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="zeroPhase">
        /// If set to <c>true</c>, the output buffer will be shifted <see cref="PhaseShift"/> samples left
        /// in order to achieve a zero phase response. Therefore, the output buffer will have 
        /// <see cref="PhaseShift"/> samples less.
        /// </param>
        public unsafe void Filter(float* input, float* output, int count, int stride, bool zeroPhase)
        {
            fixed (double* taps = _taps, sr = _sr)
            {
                int n = 0;
                if (zeroPhase)
                    for (; n < PhaseShift; n++, input += stride)
                        DoSample(taps, sr, *input);
                for (; n < count; n++, input += stride, output += stride)
                    *output = (float) DoSample(taps, sr, *input);
            }
        }

        /// <summary>
        /// Filters the specified single-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input" />.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output" />.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <exception cref="System.ArgumentOutOfRangeException">
        /// count
        /// or
        /// inputIndex
        /// or
        /// outputIndex
        /// or
        /// count - There is not enough space in the arrays
        /// or
        /// stride
        /// </exception>
        void IZeroPhaseDigitalFilter.Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count, int stride = 1)
        {
            Filter(input, inputIndex, output, outputIndex, count, stride, true);
        }

        /// <summary>
        /// Filters the specified single-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input" />.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output" />.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="zeroPhase">
        /// If set to <c>true</c>, the output buffer will be shifted <see cref="PhaseShift"/> samples left
        /// in order to achieve a zero phase response. Therefore, the output buffer will have 
        /// <see cref="PhaseShift"/> samples less.
        /// </param>
        /// <exception cref="System.ArgumentOutOfRangeException">count
        /// or
        /// inputIndex
        /// or
        /// outputIndex
        /// or
        /// count - There is not enough space in the arrays
        /// or
        /// stride</exception>
        public unsafe void Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count, int stride = 1, bool zeroPhase = false)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                zeroPhase && outputIndex + (count - PhaseShift) * stride > output.Length ||
                !zeroPhase && outputIndex + count * stride > output.Length ||
                outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
                Filter(inputPtr + inputIndex, outputPtr + outputIndex, count, stride, zeroPhase);
        }

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="zeroPhase">
        /// If set to <c>true</c>, the output buffer will be shifted <see cref="PhaseShift"/> samples left
        /// in order to achieve a zero phase response. Therefore, the output buffer will have 
        /// <see cref="PhaseShift"/> samples less.
        /// </param>
        public unsafe void Filter(double* input, double* output, int count, int stride, bool zeroPhase)
        {
            fixed (double* taps = _taps, sr = _sr)
            {
                int n = 0;
                if (zeroPhase)
                    for (; n < PhaseShift; n++, input += stride)
                        DoSample(taps, sr, *input);
                for (; n < count; n++, input += stride, output += stride)
                    *output = DoSample(taps, sr, *input);
            }
        }

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input" />.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output" />.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <exception cref="System.ArgumentOutOfRangeException">
        /// count
        /// or
        /// inputIndex
        /// or
        /// outputIndex
        /// or
        /// count - There is not enough space in the arrays
        /// or
        /// stride
        /// </exception>
        void IZeroPhaseDigitalFilter.Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count, int stride = 1)
        {
            Filter(input, inputIndex, output, outputIndex, count, stride, true);
        }

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input" />.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output" />.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="zeroPhase">
        /// If set to <c>true</c>, the output buffer will be shifted <see cref="PhaseShift"/> samples left
        /// in order to achieve a zero phase response. Therefore, the output buffer will have 
        /// <see cref="PhaseShift"/> samples less.
        /// </param>
        /// <exception cref="System.ArgumentOutOfRangeException">count
        /// or
        /// inputIndex
        /// or
        /// outputIndex
        /// or
        /// count - There is not enough space in the arrays
        /// or
        /// stride</exception>
        public unsafe void Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count, int stride = 1, bool zeroPhase = false)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                zeroPhase && outputIndex + (count - PhaseShift) * stride > output.Length || 
                !zeroPhase && outputIndex + count * stride > output.Length || 
                outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
                Filter(inputPtr + inputIndex, outputPtr + outputIndex, count, stride, zeroPhase);
        }
    }
}