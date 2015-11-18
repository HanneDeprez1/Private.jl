using EEG
using Base.Test
using Logging
using Private
using TimeModels

Logging.configure(level=DEBUG)


## Check functions work on real data

fname = joinpath(dirname(@__FILE__), "../data", "test_Hz19.5-testing.bdf")

s = read_SSR(fname)
println(s)

s.modulationrate = assr_frequency(40)

s = rereference(s, "Cz")

s = kalman_filter(s)

println(s.processing["statistics"])


## Check values using simulated data


model = Private.acoustic_model(assr_frequency(40), 8192.0, 1e-10, 5.0, 0.01, [1.0, 1.0])

new_data = Array(Float64, 8192*10, 6)
for n in 1:size(s.data, 2)
    x, new_data[:, n] = TimeModels.simulate(model, 8192*10)
end
s.data = new_data

s = kalman_filter(s)
println(s.processing["statistics"])

@test_approx_eq_eps mean(s.processing["statistics"][:SignalAmplitude][7:end]) sqrt(2) 0.2


