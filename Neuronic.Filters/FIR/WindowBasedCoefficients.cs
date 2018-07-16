using System.Collections.Generic;

namespace Neuronic.Filters.FIR
{
    public abstract class WindowBasedCoefficients : FiniteImpulseResponseCoefficients
    {
        public WindowBasedCoefficients(int n, double fs) : base(n, fs)
        {
        }

        public IWindow Window { get; set; } = FIR.Window.Hamming;

        protected void TapperEdges(IList<double> coeffs)
        {
            Window?.ApplyTo(coeffs);
        }
    }
}