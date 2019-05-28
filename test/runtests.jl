using Test
using DelimitedFiles
using SpectralEnvelope
using Random

cd(@__DIR__)
test = readdlm("DNA_data.txt")
x,y,e = spectral_envelope(test;m=0)
@test round(spectral_envelope(test)[2][5]; digits = 5) == round(0.17470002350548297;digits = 5)
@test get_mappings(0.33,x,y,e,["a","c","g","t"])[1] == "a : -0.3959646304004004"
