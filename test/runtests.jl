using Test
using DelimitedFiles
push!(LOAD_PATH,"C:\\Users\\Joy\\.julia\\dev")
using SpectralEnvelope

cd(@__DIR__)
test = readdlm("DNA_data.txt")
@test spectral_envelope(test)[2][5] == 0.17470002350548378
