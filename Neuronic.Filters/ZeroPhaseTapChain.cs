using System.Collections.Generic;

namespace Neuronic.Filters
{
    /// <summary>
    /// Represents a zero-phase finite impulse response (FIR) digital filter.
    /// </summary>
    public class ZeroPhaseTapChain : TapChain, IZeroPhaseDigitalFilter
    {
        private int _offset;

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroPhaseTapChain"/> class.
        /// </summary>
        /// <param name="taps"></param>
        public ZeroPhaseTapChain(IList<double> taps) : base(taps)
        {
        }

        /// <summary>
        /// Gets the amount of extra samples the input buffer should have if you want to obtain output buffers with a fixed size.
        /// </summary>
        public int DroppedSamples => PhaseShift - _offset;

        /// <summary>
        /// Reset's the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        public override void Reset()
        {
            base.Reset();
            _offset = 0;
        }

        /// <inheritdoc />
        protected override unsafe int ProtectedFilter(double* input, double* output, int count, int stride,
            double* taps, double* sr)
        {
            var n = 0;
            for (; _offset < PhaseShift && n < count; n++, _offset++, input += stride)
                DoSample(taps, sr, *input);
            return base.ProtectedFilter(input, output, count - n, stride, taps, sr);
        }

        /// <inheritdoc />
        protected override unsafe int ProtectedFilter(float* input, float* output, int count, int stride, double* taps,
            double* sr)
        {
            var n = 0;
            for (; _offset < PhaseShift && n < count; n++, _offset++, input += stride)
                DoSample(taps, sr, *input);
            return base.ProtectedFilter(input, output, count - n, stride, taps, sr);
        }
    }
}