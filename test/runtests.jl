using Test
using DelimitedFiles
using SpectralEnvelope
using Random

cd(@__DIR__)
test = readdlm("DNA_data.txt")
x,y,e,c = spectral_envelope(test;m=0)
@test round(spectral_envelope(test)[2][5]; digits = 5) == round(0.2819667989638535;digits = 5)
