using EEG
using Base.Test
using Logging
using Private

Logging.configure(level=DEBUG)

fname = joinpath(dirname(@__FILE__), "../data", "test_Hz19.5-testing.bdf")

s = read_ASSR(fname)

s.modulation_frequency = 40

s = rereference(s, "Cz")

s = extract_epochs(s)

s = hotelling(s)

@test_approx_eq_eps s.processing["hotelling1"][:SignalPhase][2:end]  [-2.0242503747 1.0592243553 0.1241683641 0.4413366056 -0.6632078384] 0.01

@test_approx_eq_eps s.processing["hotelling1"][:SignalPower][2:end]  [0.0245292215 0.0139141820 0.7170930945 0.0121231951  0.0419922102] 0.001

@test_approx_eq_eps s.processing["hotelling1"][:NoisePower][2:end]   [0.1625531784 0.1624553511 0.2744051671 0.1951633702  0.2178213209] 0.001

@test_approx_eq_eps s.processing["hotelling1"][:Statistic][2:end]    [0.8264963627 0.8860469460 0.0565369725 0.9547960758 0.7870649099] 0.001

@test_approx_eq_eps s.processing["hotelling1"][:SNRdB][2:end]        10*log10([0.1508996733 0.0856492686 2.6132638177 0.0621181890  0.1927828278]) 0.001

