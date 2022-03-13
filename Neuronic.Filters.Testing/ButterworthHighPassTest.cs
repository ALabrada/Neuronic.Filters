using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Transforms;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.Filters.IIR;
using Neuronic.Filters.Testing.Properties;

namespace Neuronic.Filters.Testing
{
    [TestClass]
    public class ButterworthHighPassTest
    {        
        [TestMethod]
        public void TestHighPass08()
        {
            const int order = 8;
            const double fs = 44100d;
            const double cutoffFrequency = 500d;
            const double error = 1e-4;

            var coeff = new HighPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(Resources.HighPass08).ToList();

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 2.0);

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass12()
        {
            const int order = 12;
            const double fs = 32000d;
            const double cutoffFrequency = 10d;
            const double error = 1e-4;

            var coeff = new HighPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(Resources.HighPass12).ToList();

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 1.0);

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass16()
        {
            const int order = 16;
            const double fs = 31250d;
            const double cutoffFrequency = 100d;
            const double error = 1e-4;

            var coeff = new HighPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = Helpers.LoadCsv(Resources.HighPass16).ToList();

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 1.0);

            Assert.AreEqual(expected.Count, chain.Count);
            for (int i = 0; i < expected.Count; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }

        [TestMethod]
        public void TestHighPass2()
        {
            const int order = 2;
            const double fs = 44100d;
            const double cutoffFrequency = 2000d;
            const double error = 1e-4;

            var coeff = new HighPassButterworthCoefficients(order, fs, cutoffFrequency);
            var chain = coeff.Calculate();

            var expected = new[]
            {
                new Biquad(0.817365347198867,-1.63473069439773,0.817365347198867,1,-1.60109239418362,0.668368994611848),
            };

            Assert.AreEqual(1, chain.GetTransitionBands(fs).Count());
            foreach (var f in chain.GetTransitionBands(fs))
                Assert.AreEqual(f, cutoffFrequency, 5.0);

            Assert.AreEqual(expected.Length, chain.Count);
            for (int i = 0; i < expected.Length; i++)
                Helpers.ValidateBiquad(expected[i], chain[i], error);
        }
    }
}