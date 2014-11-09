using EEG
using Private
using Logging

Logging.configure(level=WARNING)

fname = joinpath(dirname(@__FILE__), "../test/data", "test_Hz19.5-testing.bdf")

# Run everything once to prime
a = read_SSR(fname, min_epoch_length=8380, max_epoch_length=8390)
a = highpass_filter(a)
a = rereference(a, "Cz")
a = extract_epochs(a)
a = create_sweeps(a)
a = ftest(a, [1:50])
a = hotelling(a, [1:50])


println()
println("Read SSR")
@time a = read_SSR(fname, min_epoch_length=8380, max_epoch_length=8390)

println()
println("Filter")
@time a = highpass_filter(a)

println()
println("Reference")
@time a = rereference(a, "Cz")

println()
println("Epochs")
@time a = extract_epochs(a)

println()
println("Sweeps")
@time a = create_sweeps(a)

println()
println("F test")
@time a = ftest(a, [1:50])

println()
println("Hotelling")
@time a = hotelling(a, [1:50])

@profile hotelling(a, [1:50])
Profile.print()
