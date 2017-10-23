using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.Filters
{
    public class BiquadChain : IReadOnlyList<Biquad>
    {
        private readonly Biquad[] _coeffs;
        private readonly double[] _yn;
        private readonly double[] _yn1;
        private readonly double[] _yn2;
        private double _xn1;
        private double _xn2;

        public BiquadChain(IList<Biquad> coefficients, double gain)
        {
            _coeffs = new Biquad[coefficients.Count];
            coefficients.CopyTo(_coeffs, 0);
            Gain = gain;
            _yn = new double[_coeffs.Length];
            _yn1 = new double[_coeffs.Length];
            _yn2 = new double[_coeffs.Length];
        }

        public double Gain { get; }

        private void Reset()
        {
            _xn1 = 0;
            _xn2 = 0;
            Array.Clear(_yn, 0, Count);
            Array.Clear(_yn1, 0, Count);
            Array.Clear(_yn2, 0, Count);
        }

        public unsafe void Process(float* input, float* output, int count, int stride)
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

        public unsafe void Process(float[] input, int inputIndex, float[] output, int outputIndex, int count)
        {
            if (count <= 0) throw new ArgumentOutOfRangeException(nameof(count));
            if (inputIndex < 0) throw new ArgumentOutOfRangeException(nameof(inputIndex));
            if (outputIndex < 0) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (inputIndex + count > input.Length || outputIndex + count > output.Length)
                throw new ArgumentOutOfRangeException(nameof(count), "There is not enough space in the arrays");
            fixed (float* inputPtr = input, outputPtr = output)
            {
                Process(inputPtr + inputIndex, outputPtr + outputIndex, count, 1);
                Reset();
                Process(inputPtr + inputIndex + count - 1, outputPtr + outputIndex + count - 1, count, -1);
                Reset();
            }
        }

        public IEnumerator<Biquad> GetEnumerator()
        {
            return _coeffs.Cast<Biquad>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public int Count => _coeffs.Length;

        public Biquad this[int index] => _coeffs[index];
    }
}