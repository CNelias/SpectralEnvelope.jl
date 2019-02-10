# Spectral Envelope
A julia implementation of the spectral envelope method, used in categorical data analysis.

### Note : I will write a better README soon.

for a given time series TS composed of categorical data, feed it to the function spectral_envelope :
```
x,y = spectral_envelope(TS)
```
and then you can plot the results using your prefered package :
```
using Plots
plot(x,y)
```


