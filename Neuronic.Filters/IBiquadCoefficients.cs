using System.Collections.Generic;

namespace Neuronic.Filters
{
    public interface IBiquadCoefficients
    {
        double Calculate(IList<Biquad> coeffs);
    }
}