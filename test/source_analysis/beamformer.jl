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

x = repmat(collect(1:5.0), 25)
y = repmat(vec(ones(5) * collect(1:5)'), 5)
z = vec(ones(5*5) * collect(1:5)')
t = [1.0]
H = rand(125, 3, 6)
L = Leadfield(H, x, y, z, a.channel_names)

v = beamformer_lcmv(a, b, L)
v = beamformer_lcmv(a, b, L, subspace = 0.0,  regularisation = 0.0)
v = beamformer_lcmv(a, b, L, subspace = 0.0,  regularisation = 0.001)
v = beamformer_lcmv(a, b, L, subspace = 0.9,  regularisation = 0.001)
v = beamformer_lcmv(a, b, L, subspace = 0.99, regularisation = 0.001)

show(v)

