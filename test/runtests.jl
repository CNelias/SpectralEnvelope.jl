using Test
using DelimitedFiles
using SpectralEnvelope

cd(@__DIR__)
test = readdlm("DNA_data.txt")
@test round(spectral_envelope(test)[2][5]; digits = 5) == round(0.17470002350548297;digits = 5)
