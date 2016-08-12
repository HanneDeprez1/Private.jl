using EEG
using Base.Test
using Logging
using Private

fname = joinpath(dirname(@__FILE__), "../data", "test_Hz19.5-testing.bdf")

s = read_SSR(fname)

s.modulationrate = 40

s = rereference(s, "Cz")

s = extract_epochs(s)

s = hotelling(s)

# Values checked against biopil

@test_approx_eq_eps s.processing["statistics"][:SignalPhase][2:end] [-2.0242503747 1.0592243553 0.1241683641 0.4413366056 -0.6632078384] 0.01

@test_approx_eq_eps s.processing["statistics"][:SignalAmplitude][2:end] sqrt([0.0245292215 0.0139141820 0.7170930945 0.0121231951  0.0419922102]) 0.001

@test_approx_eq_eps s.processing["statistics"][:NoiseAmplitude][2:end] sqrt([0.1625531784 0.1624553511 0.2744051671 0.1951633702  0.2178213209]) 0.001

@test_approx_eq_eps s.processing["statistics"][:Statistic][2:end] [0.8264963627 0.8860469460 0.0565369725 0.9547960758 0.7870649099] 0.001

@test_approx_eq_eps s.processing["statistics"][:SNRdB][2:end] 10*log10([0.1508996733 0.0856492686 2.6132638177 0.0621181890  0.1927828278]) 0.001

# old syntax
s = hotelling(s, 20)

# Values not checked against MATLAB, just checking they dont change
@test_approx_eq_eps s.processing["statistics"][:SNRdB][2:end] [-8.213117489375172,-10.672763757064203,4.171832757032964,-12.067811128392261,-7.149317198518387,NaN,2.6057187061929428,0.31160350083134575,-0.4669094337860864,-2.8021086381963705,-4.159205467185502] 0.001
