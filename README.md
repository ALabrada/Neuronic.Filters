# Neuronic.Filters
A collection of online digital filters for Digital Signal Processing (DSP). It provides low-pass, high-pass, band-pass and band-stop implementations of the Butterworth filters (adapted from the C++ implementations made by [Ruoho Ruotsi](https://github.com/ruohoruotsi/Butterworth-Filter-Design)). 

## Code Example

```
using Neuronic.Filters.Butterwoth;
...
var coeff = new LowPassButtersworthCoefficients(8, 44100d, 500d);
var chain = coeff.Calculate();
chain.Filter(signal, 0, signal, 0, signal.Length);
```

## Instalation

The binaries are available at [Nuget](https://www.nuget.org/packages/Neuronic.Filters/). To install it run the command `Install-Package Neuronic.Filters`.

## API Reference

The library provides several implementations of the `IBiquadCoefficients` interface. They can be used to calculate the parameters of the desired IIR filter in the form of Biquad sections. Then, a `BiquadChain` should be created in order to filter the desired sample buffers. 
