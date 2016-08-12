using EEG
using Base.Test

#######################################
#
# Test beamformer
#
#######################################


fname = joinpath(dirname(@__FILE__), "..", "data", "test_Hz19.5-testing.bdf")

a = read_SSR(fname)

b = deepcopy(a)
b.data = rand(size(b.data))

a = extract_epochs(a)
b = extract_epochs(b)

#
# Fake a leadfield
#

x = repmat(collect(-2:2.0), 25)
y = repmat(vec(ones(5) * collect(-2:2.0)'), 5)
z = vec(ones(5*5) * collect(-2:2.0)')
t = [1.0]
H = rand(125, 3, 6)
L = Leadfield(H, x, y, z, channelnames(a))

v = beamformer_lcmv(a, b, L)
v = beamformer_lcmv(a, b, L, subspace = 0.0,  regularisation = 0.0)
v = beamformer_lcmv(a, b, L, subspace = 0.0,  regularisation = 0.001)
v = beamformer_lcmv(a, b, L, subspace = 0.9,  regularisation = 0.001)
v = beamformer_lcmv(a, b, L, subspace = 0.99, regularisation = 0.001)

show(v)

cpsd = cross_spectral_density(a)
@test isa(cpsd, Array{Complex{Float64}, 2})


# Test bilateral correction
v = correct_midline(v)
@test isa(v, EEG.VolumeImage)
