| **Travis**     | **Appveyor** |
|:---------------:|:-----:|
|[![Build Status](https://travis-ci.com/johncwok/SpectralEnvelope.jl.svg?branch=master)](https://travis-ci.com/johncwok/SpectralEnvelope.jl)| [![Build status](https://ci.appveyor.com/api/projects/status/q9ets366or6204u6?svg=true)](https://ci.appveyor.com/project/johncwok/spectralenvelope-jl)|




# Spectral Envelope
A fast and easy to use julia implementation of the spectral envelope method, used in categorical data analysis.

The **spectral envelope** is a tool to study cyclic behaviors in categorical data. It is more efficient than the traditional approach of attributing a different number to each category before computing the power-spectral density.<br/>

For each frequency in the spectrum, the **spectral envelope** finds an optimal real-numbered mapping that maximizes the power-spectral density at this point. Hence the name: no matter what mapping is choosen for each category, the power-spectral density will always be bounded by the spectral envelope.

The spectral envelope was defined by David S. Stoffer in *DAVID S. STOFFER, DAVID E. TYLER, ANDREW J. MCDOUGALL, Spectral analysis for categorical time series: Scaling and the spectral envelope*.\

## Usage 
The main function is:
```spectral_envelope 
spectral_envelope(ts; m = 3)

  Input
    -ts : Array containing the time series to be analysed.
    -m : Smoothing parameter. corresponds to how many neighboring points 
        are to be involved in the smoothing (weighted average). Defaults to 3.
  Returns 
    -freq : Array containing the frequency of the power-spectrum (or spectral envelope)
    -se : Values of the spectral envelope for each frequency in 'freq'.
    -eigvec : Array containing the optimal real-valued mapping for each frequency point.
    -categories : the categories which are present in the data.
```
To use the spectral envelope, call the function ```spectral_envelope```, you can then easily plot the results and extract the mapping for a given frequency.
```Julia
f, se, mappings, categories = spectral_envelope(data; m = 4)
#plotting the results
Using Plots
plot(f, se)
```
<img src=https://user-images.githubusercontent.com/34754896/91550431-d2092600-e928-11ea-8547-7fc086d41d7d.PNG width = "600">

To get the **optimal mappings** for a given frequency, you can use the ```get_mapping(data, freq; m = 3)```. With the previous DNA example, we see a peak at 0.33. To get the corresponding mappings:
```Julia
mappings = get_mappings(data, 0.33)
>> position of peak: 0.33 strengh of peak: 0.43
print(mappings)
>> ["A : 0.58", "G : -0.59", "T : 0.68", "C : 0.0"]
```
The function scans the vincinity of the provided goal frequency and returns the mapping for the found maxima. It also prints the positions and intensity of the peak so that you may control that you actually identified the desired peak and not a nearby sub-peak.<br/>
The codons A and T have an equivalent mapping and so they have similar functions (from the point of view of the time-series).

### Installation and import 
```Julia
# installing the module
Using Pkg
Pkg.clone(“https://github.com/johncwok/SpectralEnvelope.jl.git”)
# importing the module
Using SpectralEnvelope
```

## To-do
* Implement windowing & averaging (periodogram bias correction).
* Implement bootstrap confidence intervals.
