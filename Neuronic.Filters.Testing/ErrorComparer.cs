using System;
using System.Collections.Generic;

namespace Neuronic.Filters.Testing
{
    class ErrorComparer : IComparer<double>
    {
        public ErrorComparer(double error)
        {
            Error = error;
        }

        public double Error { get; }

        public int Compare(double x, double y)
        {
            if (Math.Abs(x - y) <= Error)
                return 0;
            return x.CompareTo(y);
        }
    }
}