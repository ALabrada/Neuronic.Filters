using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.Filters.Testing
{
    class FrequencySpectrum
    {
        private readonly double[] _data;
        private readonly double[] _energy;
        private readonly double _samplingFrequency;

        public FrequencySpectrum(double[] data, double samplingFrequency)
        {
            var count = data.Length / 2;
            _data = data;
            _energy = new double[data.Length];
            Array.Copy(data, _energy, _energy.Length);
            _samplingFrequency = samplingFrequency;
            var im = new double[_energy.Length];
            FourierTransform2.FFT(_energy, im, FourierTransform.Direction.Forward);
            Helpers.CalculateEnergy(_energy, im, count);
            Maximum = _energy.Take(_energy.Length / 2).Max();
            //for (int i = 0; i < _energy.Length / 2; i++)
            //    _energy[i] /= Maximum;
        }

        public double Maximum { get; }

        public double this[double frequency]
        {
            get
            {
                var fs = _samplingFrequency;
                var step = fs / _energy.Length;
                return _energy[(int) Math.Floor(frequency / step)];
            }
        }

        public IEnumerable<double> GetPeaks(double error)
        {
            var fs = _samplingFrequency;
            var step = fs / _energy.Length;
            var count = _energy.Length / 2;
            for (int i = 1; i < count - 1; i++)
            {
                var freq = i * step;
                if (_energy[i] > _energy[i - 1] && _energy[i] > _energy[i + 1] &&
                    _energy[i] >= error * Maximum)
                {
                    yield return freq;
                }
            }
        }
    }
}