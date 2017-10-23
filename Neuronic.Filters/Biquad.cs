namespace Neuronic.Filters
{
    public struct Biquad
    {
        public Biquad(double b0, double b1, double b2, double a0, double a1, double a2)
        {
            B0 = b0 / a0;
            B1 = b1 / a0;
            B2 = b2 / a0;
            A1 = -a1 / a0;
            A2 = -a2 / a0;
        }

        public double B0 { get; }

        public double B1 { get; }

        public double B2 { get; }

        public double A0 => 1;

        public double A1 { get; }

        public double A2 { get; }

        public override string ToString()
        {
            return $"B: [{B0:F4}, {B1:F4}, {B2:F4}], A: [{A0:F4}, {A1:F4}, {A2:F4}]";
        }
    }
}