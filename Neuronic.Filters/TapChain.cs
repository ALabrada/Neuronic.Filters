using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters
{
    /// <summary>
    /// Represents a linear-phase finite impulse response (FIR) digital filter.
    /// </summary>
    public class TapChain:
#if NET40
        IEnumerable<double>
#else
        IReadOnlyList<double>
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
        public virtual void Reset()
        {
            Array.Clear(_sr, 0, Count);
        }

        /// <summary>
        /// Processes the specified sample.
        /// </summary>
        /// <param name="taps">A buffer with the filter coefficients.</param>
        /// <param name="sr">A buffer with the last n processed samples.</param>
        /// <param name="dataSample">The sample to filter.</param>
        /// <returns>The processed sample.</returns>
        protected unsafe double DoSample(double* taps, double* sr, double dataSample)
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
        /// <param name="taps">A buffer with the filter coefficients.</param>
        /// <param name="sr">A buffer with the last n processed samples.</param>
        protected virtual unsafe int ProtectedFilter(float* input, float* output, int count, int stride, double* taps,
            double* sr)
        {
            for (var n = 0; n < count; n++, input += stride, output += stride)
                *output = (float)DoSample(taps, sr, *input);
            return count;
        }

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="taps">A buffer with the filter coefficients.</param>
        /// <param name="sr">A buffer with the last n processed samples.</param>
        protected virtual unsafe int ProtectedFilter(double* input, double* output, int count, int stride, double* taps,
            double* sr)
        {
            for (int n = 0; n < count; n++, input += stride, output += stride)
                *output = DoSample(taps, sr, *input);
            return count;
        }

        /// <summary>
        /// Filters the specified single-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <returns>How many samples were actually written to <paramref name="output"/>.</returns>
        public virtual unsafe int Filter(float* input, float* output, int count, int stride)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* taps = _taps, sr = _sr)
            {
                return ProtectedFilter(input, output, count, stride, taps, sr);                
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
        public virtual unsafe int Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count, int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0 || inputIndex >= input.Length) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0 || outputIndex >= output.Length) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            var inputEnd = inputIndex + count * stride;
            var outputEnd = outputIndex + count * stride;
            if (inputEnd > input.Length || inputEnd < 0 ||
                outputEnd > output.Length || outputEnd < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            fixed (float* inputPtr = input, outputPtr = output)
                return Filter(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
        }

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input" /> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public virtual unsafe int Filter(double* input, double* output, int count, int stride)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* taps = _taps, sr = _sr)
            {
                return ProtectedFilter(input, output, count, stride, taps, sr);
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
        /// <exception cref="System.ArgumentOutOfRangeException">count
        /// or
        /// inputIndex
        /// or
        /// outputIndex
        /// or
        /// count - There is not enough space in the arrays
        /// or
        /// stride</exception>
        public virtual unsafe int Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count, int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0 || inputIndex >= input.Length) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0 || outputIndex >= output.Length) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            var inputEnd = inputIndex + count * stride;
            var outputEnd = outputIndex + count * stride;
            if (inputEnd > input.Length || inputEnd < 0 ||
                outputEnd > output.Length || outputEnd < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");            
            fixed (double* inputPtr = input, outputPtr = output)
                return Filter(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
        }
        
        /// <summary>
        /// Builds a zero-phase version of this filter that works by removing the linear-phase delay.
        /// </summary>
        /// <returns>A <see cref="ZeroPhaseTapChain"/> filter with the same properties.</returns>
        public ZeroPhaseTapChain ToZeroPhase()
        {
            return new ZeroPhaseTapChain(_taps);
        }
    }
}