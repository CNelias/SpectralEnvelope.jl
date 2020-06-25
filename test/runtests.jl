using Test
using DelimitedFiles
using SpectralEnvelope
using Random

cd(@__DIR__)
test = readdlm("DNA_data.txt")
x,y,e,c = spectral_envelope(test;m=0)
@test round(spectral_envelope(test)[2][5]; digits = 5) == round(0.17517853879301779;digits = 5)
@test get_mappings(0.33,x,y,e,c)[1] == "A : 0.702"
