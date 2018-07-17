using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    /// <summary>
    /// Abstraction of a FIR filter designer that can use a windowing function to
    /// truncate the infinite impulse response.
    /// </summary>
    /// <seealso cref="Neuronic.Filters.FIR.FiniteImpulseResponseCoefficients" />
    public abstract class WindowBasedCoefficients : FiniteImpulseResponseCoefficients
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="WindowBasedCoefficients"/> class.
        /// </summary>
        /// <param name="n">The filter order.</param>
        /// <param name="fs">The sampling frequency.</param>
        public WindowBasedCoefficients(int n, double fs) : base(n, fs)
        {
        }

        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        /// <value>
        /// The window function.
        /// </value>
        public IWindow Window { get; set; } = FIR.Window.Hamming;

        /// <summary>
        /// Applies the windowing function to the specified filter coefficients.
        /// </summary>
        /// <param name="coeffs">The coefficients.</param>
        protected void ApplyWindow(IList<double> coeffs)
        {
            Window?.ApplyTo(coeffs);
        }
    }
}