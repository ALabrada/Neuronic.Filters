namespace Neuronic.Filters
{
    /// <summary>
    /// Structure that holds the coefficients of a biquad filter (see <see ref="https://en.wikipedia.org/wiki/Digital_biquad_filter"/>).
    /// </summary>
    public struct Biquad
    {
        /// <summary>
        /// Initializes a new instance of <see cref="Biquad"/>.
        /// </summary>
        /// <param name="b0">The b<sub>0</sub> coefficient.</param>
        /// <param name="b1">The b<sub>1</sub> coefficient.</param>
        /// <param name="b2">The b<sub>2</sub> coefficient.</param>
        /// <param name="a0">The a<sub>0</sub> coefficient.</param>
        /// <param name="a1">The a<sub>1</sub> coefficient.</param>
        /// <param name="a2">The a<sub>2</sub> coefficient.</param>
        public Biquad(double b0, double b1, double b2, double a0, double a1, double a2)
        {
            B0 = b0 / a0;
            B1 = b1 / a0;
            B2 = b2 / a0;
            A1 = -a1 / a0;
            A2 = -a2 / a0;
        }

        /// <summary>
        /// The b<sub>0</sub> coefficient.
        /// </summary>
        public double B0 { get; }
        /// <summary>
        /// The b<sub>1</sub> coefficient.
        /// </summary>
        public double B1 { get; }
        /// <summary>
        /// The b<sub>2</sub> coefficient.
        /// </summary>
        public double B2 { get; }
        /// <summary>
        /// The a<sub>0</sub> coefficient.
        /// </summary>
        public double A0 => 1;
        /// <summary>
        /// The a<sub>1</sub> coefficient.
        /// </summary>
        public double A1 { get; }
        /// <summary>
        /// The a<sub>2</sub> coefficient.
        /// </summary>
        public double A2 { get; }

        /// <summary>Returns the fully qualified type name of this instance.</summary>
        /// <returns>A <see cref="T:System.String" /> containing a fully qualified type name.</returns>
        /// <filterpriority>2</filterpriority>
        public override string ToString()
        {
            return $"B: [{B0:F4}, {B1:F4}, {B2:F4}], A: [{A0:F4}, {A1:F4}, {A2:F4}]";
        }
    }
}