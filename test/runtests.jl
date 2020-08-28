using Test
using DelimitedFiles
using SpectralEnvelope
using Random

cd(@__DIR__)
test = readdlm("DNA_data.txt")
x,y,e = spectral_envelope(test; m=0)
@test round(spectral_envelope(test)[2][5]; digits = 3) == round(0.126; digits = 3)
