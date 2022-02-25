using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

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

        /// <summary>
        /// Initializes a new instance of the <see cref="BiquadChain"/> class.
        /// </summary>
        /// <param name="coefficients">The list of biquad sections.</param>
        /// <param name="gain">The overall gain of the filter.</param>
        protected BiquadChain(IList<Biquad> coefficients, double gain)
        {
            _coeffs = new Biquad[coefficients.Count];
            coefficients.CopyTo(_coeffs, 0);
            Gain = gain;
        }

        /// <summary>
        /// Gets the overall gain of the filter.
        /// </summary>
        public double Gain { get; }

        /// <summary>
        /// Reset's the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        public abstract void Reset();

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
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
                Reset();
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
                Reset();
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
            {
                FilterOnce(inputPtr + inputIndex, outputPtr + outputIndex, count, stride);
                Reset();
                FilterOnce(outputPtr + outputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -stride);
                Reset();
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
            {
                FilterOnce(inputPtr, outputPtr, count, stride);
                Reset();
                FilterOnce(outputPtr + count - 1, outputPtr + count - 1, count, -stride);
                Reset();
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
            {
                FilterOnce(inputPtr, outputPtr, count, stride);
                Reset();
                FilterOnce(outputPtr + count - 1, outputPtr + count - 1, count, -stride);
                Reset();
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

        protected Biquad[] Stages => _coeffs;

        /// <summary>Gets the element at the specified index in the read-only list.</summary>
        /// <returns>The element at the specified index in the read-only list.</returns>
        /// <param name="index">The zero-based index of the element to get. </param>
        public Biquad this[int index] => Stages[index];
    }

    public class DirectFormIBiquadChain : BiquadChain
    {
        private readonly State[] _states;
        private double _ac;

        public DirectFormIBiquadChain(IList<Biquad> coefficients, double gain, double ac = 1e-8) : base(coefficients, 1.0)
        {
            _states = new State[Count];
            _ac = ac;
            Stages.Scale(gain);
        }

        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input;
            var ac = (_ac *= -1);
            result = state[0].Process(result, stage[0], ac);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

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

    public class DirectFormIIBiquadChain : BiquadChain
    {
        private readonly State[] _states;
        private double _ac;

        public DirectFormIIBiquadChain(IList<Biquad> coefficients, double gain, double ac = 1e-8) : base(coefficients, 1.0)
        {
            _states = new State[Count];
            _ac = ac;
            Stages.Scale(gain);
        }

        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input * Gain;
            var ac = (_ac *= -1);
            result = state[0].Process(result, stage[0], ac);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

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

    public class TransposedDirectFormIBiquadChain : BiquadChain
    {
        private readonly State[] _states;

        public TransposedDirectFormIBiquadChain(IList<Biquad> coefficients, double gain) : base(coefficients, 1.0)
        {
            _states = new State[Count];
            Stages.Scale(gain);
        }

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

    public class TransposedDirectFormIIBiquadChain : BiquadChain
    {
        private readonly State[] _states;
        private double _ac;

        public TransposedDirectFormIIBiquadChain(IList<Biquad> coefficients, double gain, double ac = 0.0) : base(coefficients, 1.0)
        {
            _states = new State[Count];
            _ac = ac;
            Stages.Scale(gain);
        }

        public override void Reset()
        {
            Array.Clear(_states, 0, Count);
        }

        private unsafe double FilterSample(double input, State* state, Biquad* stage)
        {
            var result = input * Gain;
            var ac = (_ac *= -1);
            result = state[0].Process(result, stage[0], ac);
            for (int i = 1; i < Count; i++)
                result = state[i].Process(result, stage[i]);
            return result;
        }

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

            public double Process(double x, Biquad s, double epsilon = 0d)
            {
                var output = _s1_1 + s.B0 * x + epsilon;
                _s1 = _s2_1 + s.B1 * x - s.A1 * output;
                _s2 = s.B2 * x - s.A2 * output;
                _s1_1 = _s1;
                _s2_1 = _s2;

                return output;
            }
        }
    }
}