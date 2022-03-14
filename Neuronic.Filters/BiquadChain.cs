using Neuronic.Filters.IIR;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Neuronic.Filters
{
    /// <summary>
    /// Uses a sequence of biquad filters to implement a higher-order filter.
    /// </summary>
    public abstract class BiquadChain : IZeroPhaseDigitalFilter
#if NET40
        , IEnumerable<Biquad>
#else
        , IReadOnlyList<Biquad>
#endif
    {
        private readonly Biquad[] _coeffs;
        private readonly double[] _transients;

        /// <summary>
        /// Initializes a new instance of the <see cref="BiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        /// <param name="padLen">The padding length. The default value is six times the amount of biquad sections.</param>
        protected BiquadChain(IList<Biquad> coefficients, double gain, int? padLen = null)
        {
            if (padLen.HasValue && padLen.Value < 0)
                throw new ArgumentOutOfRangeException(nameof(padLen));
            PaddingLength = padLen ?? coefficients.Count * 6;
            _transients = new double[PaddingLength * 2];
            _coeffs = new Biquad[coefficients.Count];
            coefficients.CopyTo(_coeffs, 0);
            Gain = gain;
        }

        /// <summary>
        /// Gets the overall gain of the filter.
        /// </summary>
        public double Gain { get; }

        /// <summary>
        /// Gets the padding length for zero phase filtering.
        /// </summary>
        public int PaddingLength { get; }

        /// <summary>
        /// Resets the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        public abstract void Reset();

        /// <summary>
        /// Resets the filter's state with a variable DC level.
        /// </summary>
        /// <param name="dc">The DC level.</param>
        public virtual void Reset(double dc)
        {
            Reset();
        }

        private unsafe double EstimateDC(double* ptr, int count)
        {
            if (count == 0)
                throw new ArgumentOutOfRangeException(nameof(count));

            var sum = 0.0;
            for (int i = count - 1; i >= 0; i--)
                sum += ptr[i];
            return sum / count;
        }

        /// <summary>
        /// Reads the samples used as signal padding.
        /// </summary>
        /// <param name="inputPtr">The input pointer.</param>
        /// <param name="count">The signal length.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="outputPtr">The output pointer.</param>
        /// <returns>The amount of padded samples on each side.</returns>
        protected virtual unsafe int ReadTransients(float* inputPtr, int count, int stride, double* outputPtr)
        {
            var nfact = Math.Min(_transients.Length / 2, count);
            if (nfact == 0) return 0;
            var current = outputPtr;
            var f = 2 * inputPtr[0];
            for (int i = nfact; i > 0; i--, current++)
                *current = f - inputPtr[i];
            f = 2 * inputPtr[count - 1];
            for (int i = count - nfact - 1; i < count - 1; i++, current++)
                *current = f - inputPtr[i];
            return nfact;
        }

        /// <summary>
        /// Reads the samples used as signal padding.
        /// </summary>
        /// <param name="inputPtr">The input pointer.</param>
        /// <param name="count">The signal length.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <param name="outputPtr">The output pointer.</param>
        /// <returns>The amount of padded samples on each side.</returns>
        protected virtual unsafe int ReadTransients(double* inputPtr, int count, int stride, double* outputPtr)
        {
            var nfact = Math.Min(_transients.Length / 2, count);
            if (nfact == 0) return 0;
            var current = outputPtr;
            var f = 2 * inputPtr[0];
            for (int i = nfact; i > 0; i--, current++)
                *current = f - inputPtr[i * stride];
            f = 2 * inputPtr[count - 1];
            for (int i = count - nfact - 1; i < count - 1; i++, current++)
                *current = f - inputPtr[i * stride];
            return nfact;
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public abstract unsafe void FilterOnce(float* input, float* output, int count, int stride);

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
        public unsafe int Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count,
            int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
            fixed (double* tranPtr = _transients)
            {
                var nfact = ReadTransients(inputPtr, count, stride, tranPtr);

                // Init state
                Reset(nfact <= 0 ? inputPtr[0] : EstimateDC(tranPtr, nfact));
                FilterOnce(tranPtr, tranPtr, nfact, 1);
                // Filter forward
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);

                // Init state
                Reset(nfact <= 0 ? outputPtr[outputIndex + count - 1] : EstimateDC(tranPtr + nfact, nfact));
                FilterOnce(tranPtr + nfact, tranPtr + nfact, nfact, 1);
                // Filter back
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
            }

            return count;
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
        public abstract unsafe void FilterOnce(double* input, double* output, int count, int stride);

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
        public unsafe int Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count,
            int stride = 1)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count * stride > input.Length || inputIndex + count * stride < 0 ||
                outputIndex + count * stride > output.Length || outputIndex + count * stride < 0)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
            fixed (double* tranPtr = _transients)
            {
                var nfact = ReadTransients(inputPtr, count, stride, tranPtr);

                // Init state
                Reset(nfact <= 0 ? inputPtr[0] : EstimateDC(tranPtr, nfact));
                FilterOnce(tranPtr, tranPtr, nfact, 1);
                // Filter forward
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);

                // Init state
                Reset(nfact <= 0 ? outputPtr[outputIndex + count - 1] : EstimateDC(tranPtr + nfact, nfact));
                FilterOnce(tranPtr + nfact, tranPtr + nfact, nfact, 1);
                // Filter back
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
            }

            return count;
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

#if !NET40
        /// <summary>
        /// Executes two sweeps of the filter (forward and backward). This is a zero-phase filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <remarks>
        /// This method will reset the filter.
        /// </remarks>
        public unsafe int Filter(ReadOnlySpan<float> input, Span<float> output, int count, int stride)
        {
            if (count <= 0 || input.Length < count || output.Length < count) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
            fixed (double* tranPtr = _transients)
            {
                var nfact = ReadTransients(inputPtr, count, stride, tranPtr);

                // Init state
                Reset();
                FilterOnce(tranPtr, tranPtr, nfact, 1);
                // Filter forward
                FilterOnce(inputPtr, outputPtr, count, stride);

                // Init state
                Reset();
                FilterOnce(tranPtr, tranPtr, _transients.Length, 1);
                // Filter back
                FilterOnce(outputPtr + count - 1, outputPtr + count - 1, count, -stride);
            }

            return count;
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(ReadOnlySpan<float> input, Span<float> output, int count,
            int stride = 1)
        {
            if (count <= 0 || input.Length < count || output.Length < count) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (float* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr, outputPtr, count, stride);
            }
        }

        /// <summary>
        /// Executes two sweeps of the filter (forward and backward). This is a zero-phase filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        /// <remarks>
        /// This method will reset the filter.
        /// </remarks>
        public unsafe int Filter(ReadOnlySpan<double> input, Span<double> output, int count,
            int stride = 1)
        {
            if (count <= 0 || input.Length < count || output.Length < count) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
            fixed (double* tranPtr = _transients)
            {
                var nfact = ReadTransients(inputPtr, count, stride, tranPtr);

                // Init state
                Reset();
                FilterOnce(tranPtr, tranPtr, nfact, 1);
                // Filter forward
                FilterOnce(inputPtr, outputPtr, count, stride);

                // Init state
                Reset();
                FilterOnce(tranPtr, tranPtr, _transients.Length, 1);
                // Filter back
                FilterOnce(outputPtr + count - 1, outputPtr + count - 1, count, -stride);
            }

            return count;
        }

        /// <summary>
        /// Executes one sweep of the filter.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        public unsafe void FilterOnce(ReadOnlySpan<double> input, Span<double> output, int count,
            int stride = 1)
        {
            if (count <= 0 || input.Length < count || output.Length < count) throw new ArgumentOutOfRangeException(nameof(count));
            if (stride == 0) throw new ArgumentOutOfRangeException(nameof(stride));
            fixed (double* inputPtr = input, outputPtr = output)
            {
                FilterOnce(inputPtr, outputPtr, count, stride);
            }
        }
#endif

        /// <summary>Returns an enumerator that iterates through the collection.</summary>
        /// <returns>An enumerator that can be used to iterate through the collection.</returns>
        /// <filterpriority>1</filterpriority>
        public IEnumerator<Biquad> GetEnumerator()
        {
            return Stages.Cast<Biquad>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        /// <summary>Gets the number of elements in the collection.</summary>
        /// <returns>The number of elements in the collection. </returns>
        public int Count => Stages.Length;

        /// <summary>
        /// Gets the filter's second order stages.
        /// </summary>
        protected Biquad[] Stages => _coeffs;

        /// <summary>Gets the element at the specified index in the read-only list.</summary>
        /// <returns>The element at the specified index in the read-only list.</returns>
        /// <param name="index">The zero-based index of the element to get. </param>
        public Biquad this[int index] => Stages[index];

        /// <summary>
        /// Gets the impulse response of the filter at the specified frequency.
        /// </summary>
        /// <param name="normalizedFrequency">The normalized frequency.</param>
        /// <returns>The impulse response.</returns>
        public Complex GetResponse(double normalizedFrequency)
        {
            return Stages.GetResponse(normalizedFrequency);
        }
    }

    /// <summary>
    /// A biquad filter in Direct Form I.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.BiquadChain" />
    public class DirectFormIBiquadChain : BiquadChain
    {
        private readonly State[] _states;

        /// <summary>
        /// Initializes a new instance of the <see cref="DirectFormIBiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        /// <param name="padLen">The padding length. The default value is six times the amount of biquad sections.</param>
        public DirectFormIBiquadChain(IList<Biquad> coefficients, double gain, int? padLen = null) : base(coefficients, gain, padLen)
        {
            _states = new State[Count];
        }

        /// <inheritdoc />
        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input;
            result = state[0].Process(result, stage[0]);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(float* input, float* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = (float) FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(double* input, double* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        struct State
        {
            private double _x1, _x2, _y1, _y2;

            public double Process(double x, Biquad s, double epsilon = 0d)
            {
                var output = s.B0 * x
                    + s.B1 * _x1 + s.B2 * _x2
                    - s.A1 * _y1 - s.A2 * _y2 + epsilon;
                _x2 = _x1;
                _y2 = _y1;
                _x1 = x;
                _y1 = output;

                return output;
            }
        }
    }

    /// <summary>
    /// A biquad filter in Direct Form II.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.BiquadChain" />
    public class DirectFormIIBiquadChain : BiquadChain
    {
        private readonly State[] _states;

        /// <summary>
        /// Initializes a new instance of the <see cref="DirectFormIIBiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        /// <param name="padLen">The padding length. The default value is six times the amount of biquad sections.</param>
        public DirectFormIIBiquadChain(IList<Biquad> coefficients, double gain, int? padLen = null) : base(coefficients, gain, padLen)
        {
            _states = new State[Count];
        }

        /// <inheritdoc />
        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input * Gain;
            result = state[0].Process(result, stage[0]);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(float* input, float* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = (float)FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(double* input, double* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        struct State
        {
            private double _v1, _v2;

            public double Process(double x, Biquad s, double epsilon = 0d)
            {
                var w = x - s.A1 * _v1 - s.A2 * _v2 + epsilon;
                var output = s.B0 * w + s.B1 * _v1 + s.B2 * _v2;

                _v2 = _v1;
                _v1 = w;

                return output;
            }
        }
    }

    /// <summary>
    /// A biquad filter in Transposed Direct Form I.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.BiquadChain" />
    public class TransposedDirectFormIBiquadChain : BiquadChain
    {
        private readonly State[] _states;

        /// <summary>
        /// Initializes a new instance of the <see cref="TransposedDirectFormIBiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        /// <param name="padLen">The padding length. The default value is six times the amount of biquad sections.</param>
        public TransposedDirectFormIBiquadChain(IList<Biquad> coefficients, double gain, int? padLen = null) : base(coefficients, gain, padLen)
        {
            _states = new State[Count];
        }

        /// <inheritdoc />
        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input * Gain;
            result = state[0].Process(result, stage[0]);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(float* input, float* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = (float)FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(double* input, double* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        struct State
        {
            private double _v, _s1, _s1_1, _s2, _s2_1, _s3, _s3_1, _s4, _s4_1;

            public double Process(double x, Biquad s)
            {
                // can be: in += m_s1_1;
                _v = x +_s1_1;
                var output = s.B0 * _v + _s3_1;
                _s1 = _s2_1 - s.A1 * _v;
                _s2 = -s.A2 * _v;
                _s3 = s.B1 * _v + _s4_1;
                _s4 = s.B2 * _v;

                _s4_1 = _s4;
                _s3_1 = _s3;
                _s2_1 = _s2;
                _s1_1 = _s1;

                return output;
            }
        }
    }

    /// <summary>
    /// A biquad filter in Transposed Direct Form II.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.BiquadChain" />
    public class TransposedDirectFormIIBiquadChain : BiquadChain
    {
        private readonly State[] _states;
        private bool _estimateDC;

        /// <summary>
        /// Initializes a new instance of the <see cref="TransposedDirectFormIIBiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        /// <param name="padLen">The padding length. The default value is six times the amount of biquad sections.</param>
        /// <param name="estimateDC">Whether to estimate the DC level before filtering and to initialize the state accordingly.</param>
        public TransposedDirectFormIIBiquadChain(IList<Biquad> coefficients, double gain, int? padLen = null, bool estimateDC = false) : base(coefficients, gain, padLen)
        {
            _estimateDC = estimateDC;
            _states = new State[Count];
        }

        /// <inheritdoc />
        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        /// <inheritdoc />
        public override void Reset(double dc)
        {
            if (!_estimateDC || dc.Equals(0d))
            {
                Reset();
                return;
            }

            var scale = 1.0;
            for (int i = 0; i < Stages.Length; i++)
            {
                var s = Stages[i];
                _states[i] = State.FromDCLevel(s, dc) * scale;
                scale *= (s.B0 + s.B1 + s.B2) / (s.A0 + s.A1 + s.A2);
            }
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input * Gain;
            result = state[0].Process(result, stage[0]);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(float* input, float* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = (float)FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        /// <inheritdoc />
        public override unsafe void FilterOnce(double* input, double* output, int count, int stride)
        {
            fixed (Biquad* stage = Stages)
            fixed (State* state = _states)
            {
                for (int n = 0; n < count; n++)
                {
                    *output = FilterSample(*input, state, stage);

                    input += stride;
                    output += stride;
                }
            }
        }

        struct State
        {
            private double _s1, _s1_1, _s2, _s2_1;

            private State(double z0, double z1)
            {
                _s1 = _s1_1 = z0;
                _s2 = _s2_1 = z1;
            }

            public double Process(double x, Biquad s, double epsilon = 0d)
            {
                var output = _s1_1 + s.B0 * x + epsilon;
                _s1 = _s2_1 + s.B1 * x - s.A1 * output;
                _s2 = s.B2 * x - s.A2 * output;
                _s1_1 = _s1;
                _s2_1 = _s2;

                return output;
            }

            public static State FromDCLevel(Biquad stage, double dc)
            {
                var idMinA0 = 1 + stage.A1;
                var idMinA1 = stage.A2;

                var b1 = stage.B1 - stage.A1 * stage.B0;
                var b2 = stage.B2 - stage.A2 * stage.B0;

                var aSum = 1.0;
                var cSum = 0.0;
                var z0 = (b1 + b2) / (idMinA0 + idMinA1);

                aSum += stage.A1;
                cSum += b1;
                var z1 = aSum * z0 - cSum;

                return new State(z0 * dc, z1 * dc);
            }

            public static State operator *(State state, double scale)
            {
                return new State(state._s1_1 * scale, state._s2_1 * scale);
            }
        }
    }
}