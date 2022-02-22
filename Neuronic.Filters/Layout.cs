using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace Neuronic.Filters
{
    class Layout : IList<PoleZeroPair>
    {
        private readonly List<PoleZeroPair> _analogPairs;

        public Layout(List<PoleZeroPair> analogPairs)
        {
            _analogPairs = analogPairs ?? throw new ArgumentNullException(nameof(analogPairs));
        }

        public Layout(int size) : this (new List<PoleZeroPair>(size))
        {

        }

        public PoleZeroPair this[int index] { get => ((IList<PoleZeroPair>)_analogPairs)[index]; set => ((IList<PoleZeroPair>)_analogPairs)[index] = value; }

        public int Count => ((ICollection<PoleZeroPair>)_analogPairs).Count;

        public int NumberOfPoles { get; private set; }

        public bool IsReadOnly => ((ICollection<PoleZeroPair>)_analogPairs).IsReadOnly;

        public double NormalW { get; set; }

        public double NormalGain { get; set; }

        void ICollection<PoleZeroPair>.Add(PoleZeroPair item)
        {
            ((ICollection<PoleZeroPair>)_analogPairs).Add(item);
            NumberOfPoles += 2;
        }

        public void AddPoleZeroConjugatePairs(Complex pole, Complex zero)
        {
            _analogPairs.Add(new PoleZeroPair(pole, zero, Complex.Conjugate(pole), Complex.Conjugate(zero)));
            NumberOfPoles += 2;
        }

        public void Add(Complex pole, Complex zero)
        {
            _analogPairs.Add(new PoleZeroPair(pole, zero));
            NumberOfPoles++;
        }

        public void Add(ComplexPair poles, ComplexPair zeros)
        {
            _analogPairs.Add(new PoleZeroPair(poles.First, zeros.First, poles.Second, zeros.Second));
            NumberOfPoles += 2;
        }

        public void Clear()
        {
            ((ICollection<PoleZeroPair>)_analogPairs).Clear();
        }

        public bool Contains(PoleZeroPair item)
        {
            return ((ICollection<PoleZeroPair>)_analogPairs).Contains(item);
        }

        public void CopyTo(PoleZeroPair[] array, int arrayIndex)
        {
            ((ICollection<PoleZeroPair>)_analogPairs).CopyTo(array, arrayIndex);
        }

        public IEnumerator<PoleZeroPair> GetEnumerator()
        {
            return ((IEnumerable<PoleZeroPair>)_analogPairs).GetEnumerator();
        }

        public int IndexOf(PoleZeroPair item)
        {
            return ((IList<PoleZeroPair>)_analogPairs).IndexOf(item);
        }

        public void Insert(int index, PoleZeroPair item)
        {
            ((IList<PoleZeroPair>)_analogPairs).Insert(index, item);
        }

        public bool Remove(PoleZeroPair item)
        {
            return ((ICollection<PoleZeroPair>)_analogPairs).Remove(item);
        }

        public void RemoveAt(int index)
        {
            ((IList<PoleZeroPair>)_analogPairs).RemoveAt(index);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return ((IEnumerable)_analogPairs).GetEnumerator();
        }
    }

    struct PoleZeroPair
    {
        public PoleZeroPair(Complex pole, Complex zero) : this(pole, zero, Complex.Zero, Complex.Zero)
        {
        }

        public PoleZeroPair(
            Complex pole1, Complex zero1,
            Complex pole2, Complex zero2)
        {
            Poles = new ComplexPair(pole1, pole2);
            Zeros = new ComplexPair(zero1, zero2);
        }

        public ComplexPair Poles { get; }

        public ComplexPair Zeros { get; }

        public bool IsSignlePole => Poles.Second.Equals(Complex.Zero) && Zeros.Second.Equals(Complex.Zero);
    }

    struct ComplexPair
    {
        public ComplexPair(Complex first, Complex second)
        {
            First = first;
            Second = second;
        }

        public Complex First { get; set; }

        public Complex Second { get; set; }
    }
}
