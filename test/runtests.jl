using Test
using DelimitedFiles
using SpectralEnvelope
using Random

cd(@__DIR__)
test = readdlm("DNA_data.txt",',')[1,:]
x,y,e = spectral_envelope(test; m=0)
@test round(spectral_envelope(test)[2][5]; digits = 3) == round(0.003; digits = 3)
@test round(get_mappings(test, 0.33)["A"]; digits = 2) == round(0.54; digits = 2)
