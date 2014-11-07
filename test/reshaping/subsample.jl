using EEG
using Base.Test
using Logging
using Private
Logging.configure(level=DEBUG)

fname = joinpath(dirname(@__FILE__), "../data", "test_Hz19.5-testing.bdf")

a = read_SSR(fname)

a.processing["Carrier_Frequency"] = 500

a = subsample(a, plot=false, plot_channel=4)

a = highpass_filter(a)

a = rereference(a, "Cz")

a = extract_epochs(a)

a = hotelling(a)

println(a.processing["hotelling1"])

# These values are not validated, but just to ensure it doesn't change
# TODO get true values to check with
@test_approx_eq_eps a.processing["hotelling1"][:SNRdB] [NaN, -0.0474418, -1.42819, -0.334432, -2.56155, -1.39146] 0.002
