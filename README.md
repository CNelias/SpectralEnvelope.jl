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
<img src=https://user-images.githubusercontent.com/34754896/81937423-e5031f00-95f3-11ea-986d-bb5a3689639f.png width = "600">

To get the **optimal mappings** for a given frequency more easily, you can use the ```get_mapping(goal, f, se, mappings, categories)``` function (you can use spalting for a more concise call) :
```Julia
f,se,mappings,categories =spectral_envelope(data; m =0)
get_mapping(0.33,f,se,mappings,categories)
# using spalting : get_mapping(goal, spectral_envelope(data;m=0)…)
position of peak: 0.33 strengh of peak: 0.86
["A : 0.702", "G : 0.702", "T : 0.561", "C : 0.233"]
```
The function scans the vincinity of the provided goal frequency and returns the mapping for the found maxima. It also prints the positions and intensity of the peak so that you may control that you actually identified the desired peak and not a nearby sub-peak.<br/>
In this example, we see that at the frequency ~0.33, the codons A and G have an equivalent mapping and so they have the same function (from the point of view of the time-series).

### Installation and import 
```Julia
# installing the module
Using Pkg
Pkg.clone(“https://github.com/johncwok/SpectralEnvelope.jl.git”)
# importing the module
Using SpectralEnvelope
```

## To-do
* Implement detrending and windowing for bias correction.
* Make the extraction of optimal mappings more user-friendly
