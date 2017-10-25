using System.Collections.Generic;

namespace Neuronic.Filters
{
    /// <summary>
    /// Abstraction of a class that calculates filtering coefficients for a <see cref="BiquadChain"/>.
    /// </summary>
    public interface IBiquadCoefficients
    {
        /// <summary>
        /// Calculates a set of coefficients for a filter.
        /// </summary>
        /// <param name="coeffs">The container for the coefficients.</param>
        /// <returns>The gain of the filter.</returns>
        double Calculate(IList<Biquad> coeffs);
    }
}