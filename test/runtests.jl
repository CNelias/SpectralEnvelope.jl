using Test
using DelimitedFiles
using SpectralEnvelope

cd(@__DIR__)
test = readdlm("DNA_data.txt")
@test signif(spectral_envelope(test)[2][5],4) == signif(0.17470002350548297,4)
