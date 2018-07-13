namespace Neuronic.Filters
{
    /// <summary>
    /// Abstraction of a zero-phase digital filter.
    /// </summary>
    public interface IZeroPhaseDigitalFilter
    {
        /// <summary>
        /// Reset's the filter's state. Use this if the next buffer that will be processed is not continuous after the last one.
        /// </summary>
        void Reset();

        /// <summary>
        /// Filters the specified single-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        void Filter(float[] input, int inputIndex, float[] output, int outputIndex, int count, int stride = 1);

        /// <summary>
        /// Filters the specified double-precision sample buffer.
        /// </summary>
        /// <param name="input">The source buffer.</param>
        /// <param name="inputIndex">The starting index in <paramref name="input"/>.</param>
        /// <param name="output">The destination buffer. (Use the same value of <paramref name="input"/> to execute in place).</param>
        /// <param name="outputIndex">The starting index in <paramref name="output"/>.</param>
        /// <param name="count">The number of samples to process.</param>
        /// <param name="stride">The length of the processing step. Can be used to bypass samples.</param>
        void Filter(double[] input, int inputIndex, double[] output, int outputIndex, int count, int stride = 1);
    }
}