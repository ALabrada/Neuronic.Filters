using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;

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

    public abstract class LeastSquareCoefficients : WindowBasedCoefficients
    {
        private readonly double[] _frequencies;
        private readonly double[] _amplitudes;

        public LeastSquareCoefficients(int n, double fs, IList<double> ff, IList<double> aa) : base(n, fs)
        {
            _frequencies = new double[ff.Count];
            ff.CopyTo(_frequencies, 0);
            Array.Sort(_frequencies);
            _amplitudes = new double[aa.Count];
            aa.CopyTo(_amplitudes, 0);
        }

        public bool UseScaling { get; set; } = true;

        public static double NormalizeFrequency(double f, double fs)
        {
            return 0.5 * f / fs;
        }

        private static double[] mldivide(double[,] a, IList<double> b)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Linear-phase FIR filter design using least-squares error minimization.
        /// </summary>
        /// <param name="n">The desired filter order.</param>
        /// <param name="frequencies">
        /// vector of frequency band edges in pairs, in ascending
        /// order between 0 and 1. 1 corresponds to the Nyquist frequency or half
        /// the sampling frequency.
        /// </param>
        /// <param name="amplitudes">
        /// A is a real vector the same size as F
        /// which specifies the desired amplitude of the frequency response of the
        /// resultant filter B.
        /// </param>
        /// <param name="weights">
        /// has one entry per band (so it is half the length of F and A) which
        /// tells <see cref="FIRLS"/> how much emphasis to put on minimizing the integral squared error
        /// in each band relative to the other bands.
        /// </param>
        /// <returns>
        /// a length <paramref name="n"/>+1 linear phase (real, symmetric
        /// coefficients) FIR filter which has the best approximation to the
        /// desired frequency response described by <paramref name="frequencies"/> 
        /// and <paramref name="amplitudes"/> in the least squares sense.
        /// </returns>
        internal static double[] FIRLS(int n, IList<double> frequencies, IList<double> amplitudes, IList<double> weights = null)
        {
            var pi = Math.PI;

            double sqr(double x) => x * x;
            double sinc(double t) => t == 0 ? 1 : Math.Sin(pi * t) / (pi * t);

            var vF = frequencies.ToList();
            var vM = amplitudes;
            var vW = weights;

            Debug.Assert(vF.All(x => x >= 0 && x <= 1), "The frequencies must be in [0,1].");
            Debug.Assert(vF.Count % 2 == 0, "The frequency vector must have even length.");
            Debug.Assert(vF.Count == vM.Count, "The frequency and amplitude vectors must have equal length.");

            if (vW == null)
            {
                // W = ones(length(F)/2,1);
                var newW = new List<double>(vF.Count / 2);
                newW.AddRange(Enumerable.Repeat(1d, vF.Count / 2));
                vW = newW;
            }

            // Check for valid filter length
            if (vM[vM.Count - 1] != 0 && vF[vF.Count - 1] == 1 && n % 2 == 1)
                n++;

            // N = N+1;                   % filter length
            n++;
            // F=F(:)/2;  M=M(:);  W=sqrt(W(:));  % make these guys columns
            for (int i = 0; i < vF.Count; i++)
                vF[i] /= 2d;
            //for (int i = 0; i < w.Count; i++)
            //    w[i] = Math.Sqrt(w[i]);

            var dF = vF.Zip(vF.Skip(1), (pred, suc) => suc - pred).ToList();
            Debug.Assert(dF.Any(x => x <= 0), "The frequency vector must be in order.");

            /*
            if all(dF(2:2:length(dF)-1)==0) && length(dF) > 1,
                fullband = 1;
            else
                fullband = 0;
            end
            */
            var fullband = dF.Count > 1 && Enumerable.Range(0, dF.Count).Select(x => 2 * x + 1)
                               .TakeWhile(x => x < dF.Count - 1).All(x => dF[x] == 0);

            // L=(N-1)/2;
            var l = (n - 1) / 2;

            // Nodd = rem(N,2);
            var nOdd = n % 2 != 0;

            /*
             if ~Nodd
                m=(0:L)+.5;   % type II
            else
                m=(0:L);      % type I
            end
             */
            var firstK = nOdd ? 0.0 : 0.5;
            // need_matrix = (~fullband) || (~constant_weights);
            var needMatrix = !fullband;
            /*
            if need_matrix
                I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
                I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
                G=zeros(size(I1));
            end             
             */

            var g = new double[l + 1, l + 1];
            IList<double> b = new double[l + 1];

            // for s=1:2:length(F),
            for (int s = 0; s < vF.Count; s += 2)
            {
                // m=(M(s+1)-M(s))/(F(s+1)-F(s));    %  slope
                var m = (vM[s + 1] - vM[s]) / (vF[s + 1] - vF[s]);
                // b1=M(s)-m*F(s);                   %  y-intercept
                var b1 = vM[s] - m * vF[s];
                var abs = Math.Abs(sqr(vW[(s + 1) / 2]));
                if (nOdd)
                    // b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)))...
                    //      * abs(W((s + 1) / 2) ^ 2);
                    b[0] += (b1 * (vF[s + 1] - vF[s]) + m / 2 * (vF[s + 1] * vF[s + 1] - vF[s] * vF[s]))
                            * abs;
                /*
                b = b+(m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k))...
                    * abs(W((s+1)/2)^2);
                b = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
                    - F(s)*(m*F(s)+b1)*sinc(2*k*F(s))) ...
                    * abs(W((s+1)/2)^2);
                 */
                
                for (int i = nOdd ? 1 : 0; i < b.Count; i++)
                {
                    var k = firstK + i;
                    b[i] += (m / (4 * pi * pi) * (Math.Cos(2 * pi * k * vF[s + 1]) - Math.Cos(2 * pi * k * vF[s])) / sqr(k))
                            * abs;
                    b[i] += (vF[s + 1] * (m * vF[s + 1] + b1) * sinc(2 * k * vF[s + 1])
                             - vF[s] * (m * vF[s] + b1) * sinc(2 * k * vF[s]))
                            * abs;
                }
                if (needMatrix)
                {
                    /*
                    G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))+sinc(2*I2*F(s+1))) ...
                        - .5*F(s)*(sinc(2*I1*F(s))+sinc(2*I2*F(s))) ) ...
                        * abs(W((s+1)/2)^2);
                    */
                    for (int i = 0; i <= l; i++)
                    for (int j = 0; j <= l; j++)
                    {
                        var mI1 = i + j;
                        var mI2 = i - j;
                        if (nOdd)
                            mI1++;
                        g[i, j] += (.5 * vF[s + 1] * (sinc(2 * mI1 * vF[s + 1]) + sinc(2 * mI2 * vF[s + 1]))
                                    - .5 * vF[s] * (sinc(2 * mI1 * vF[s]) + sinc(2 * mI2 * vF[s])))
                                   * abs;
                    }
                }
            }
            IList<double> a;
            if (needMatrix)
                // a=G\b;
                a = mldivide(g, b);
            else
            {
                a = new double[b.Count];
                // a=(W(1)^2)*4*b;
                for (int i = 0; i < b.Count; i++)
                    a[i] = sqr(vW[0]) * 4 * b[i];
                if (nOdd)
                    // a(1) = a(1)/2;
                    a[0] /= 2;
            }
            var h = new List<double>();
            if (nOdd)
            {
                // h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
                for (int i = l; i >= 1; i--)
                    h.Add(0.5 * a[i]);
                h.Add(a[0]);
                for (int i = 1; i <= l; i++)
                    h.Add(0.5 * a[i]);
            }
            else
            {
                // h=.5*[flipud(a); a].';
                for (int i = l; i >= 0; i--)
                    h.Add(0.5 * a[i]);
                for (int i = 0; i <= l; i++)
                    h.Add(0.5 * a[i]);
            }
            return h.ToArray();
        }

        public override void Calculate(IList<double> coeffs)
        {
            /*
            % Parse optional input arguments
            [Ftype,Wind,SCALING] = parseoptargs(Wn,varargin{:});

            % Compute the frequency vector
            [nbands,ff,Ftype] = desiredfreq(Wn,Ftype);

            % Compute the magnitude vector
            [aa,First_Band] = desiredmag(Ftype,nbands);

            % Check for appropriate filter order, increase when necessary
            [N,msg1,msg2,msgobj] = firchk(N,ff(end),aa,hilbert);
            if ~isempty(msg1), error(msgobj); end
            if ~isempty(msg2), warning(msgobj); end

            % Work with filter length (= order + 1)
            L = N + 1;

            % Check for valid window, or assign default if empty
            Wind = chkwindow(Wind,L);

            % Compute unwindowed impulse response
            if hilbert
                hh = firls(L-1,ff,aa,'h');
            else
                hh = firls(L-1,ff,aa);
            end

            % Window impulse response to get the filter
            b = hh.*Wind(:)';
            a = 1;

            if SCALING,
                % Scale so that passband is approx 1
                b = scalefilter(b,First_Band,ff,L);
            end
            */

            var hh = FIRLS(FilterOrder, _frequencies, _amplitudes);
            foreach (var x in hh)
                coeffs.Add(x);

            TapperEdges(coeffs);

            if (UseScaling)
                ScaleFilter(coeffs);
        }

        /// <summary>
        /// Scale fitler to have passband approx. equal to one.
        /// </summary>
        protected virtual void ScaleFilter(IList<double> b)
        {
            var firstBand = _amplitudes[0] > 0;
            var ff = _frequencies;
            var length = FilterOrder + 1;

            double scaleFactor;
            if (firstBand)
            {
                // b = b / sum(b);  % unity gain at DC
                scaleFactor = b.Sum();
            }
            else
            {
                /*
                if ff(4)==1
                    % unity gain at Fs/2
                    f0 = 1;
                else
                    % unity gain at center of first passband
                    f0 = mean(ff(3:4));
                end
                 */
                var f0 = ff[3] == 1 ? 1 : (ff[2] + ff[3]) / 2d;

                // b = b / abs( exp(-1i*2*pi*(0:L-1)*(f0/2))*(b.') );
                Complex sum = 0;
                for (int i = 0; i < length; i++)
                {
                    var x = new Complex(0, -1) * 2 * Math.PI * i * (f0 / 2d);
                    var exp = Complex.Exp(x);
                    sum += exp * b[i];
                }
                scaleFactor = Complex.Abs(sum);
            }
            for (int i = 0; i < b.Count; i++)
                b[i] /= scaleFactor;
        }
    }

    public class LowPassLeastSquareCoefficients : LeastSquareCoefficients
    {
        public LowPassLeastSquareCoefficients(int n, double fs, double fx) 
            : base(n, fs, 
                  new [] {0d, NormalizeFrequency(fx, fs), NormalizeFrequency(fx, fs), 1d}, 
                  new [] {1d, 1d, 0d, 0d})
        {
            if (fx <= 0 || fx >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(fx));

            CutoffFrequency = fx;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }
    }

    public class HighPassLeastSquareCoefficients : LeastSquareCoefficients
    {
        public HighPassLeastSquareCoefficients(int n, double fs, double fx)
            : base(n, fs,
                new[] { 0d, NormalizeFrequency(fx, fs), NormalizeFrequency(fx, fs), 1d },
                new[] { 0d, 0d, 1d, 1d })
        {
            if (fx <= 0 || fx >= fs / 2)
                throw new ArgumentOutOfRangeException(nameof(fx));

            CutoffFrequency = fx;
        }

        /// <summary>
        /// The cut-off frequency.
        /// </summary>
        public double CutoffFrequency { get; }
    }

    public class BandPassLeastSquareCoefficients : LeastSquareCoefficients
    {
        public BandPassLeastSquareCoefficients(int n, double fs, double f1, double f2) 
            : base(n, fs, 
                  new [] {0d, NormalizeFrequency(f1, fs), NormalizeFrequency(f1, fs), NormalizeFrequency(f2, fs), NormalizeFrequency(f2, fs), 1d },
                  new [] {0d, 0d, 1d, 1d, 0d, 0d})
        {
            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }
    }

    public class BandStopLeastSquareCoefficients : LeastSquareCoefficients
    {
        public BandStopLeastSquareCoefficients(int n, double fs, double f1, double f2)
            : base(n, fs,
                new[] { 0d, NormalizeFrequency(f1, fs), NormalizeFrequency(f1, fs), NormalizeFrequency(f2, fs), NormalizeFrequency(f2, fs), 1d },
                new[] { 1d, 1d, 0d, 0d, 1d, 1d })
        {
            if (f1 > f2)
            {
                var temp = f2;
                f2 = f1;
                f1 = temp;
            }
            FirstCutoffFrequency = f1;
            SecondCutoffFrequency = f2;
        }

        /// <summary>
        /// The minor cut-off frequency.
        /// </summary>
        public double FirstCutoffFrequency { get; }
        /// <summary>
        /// The major cut-off frequency.
        /// </summary>
        public double SecondCutoffFrequency { get; }
    }
}