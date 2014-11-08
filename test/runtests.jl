using Private
using Base.Test
using Logging

Logging.configure(level=DEBUG)

include("statistics/hotelling.jl")
include("reshaping/subsample.jl")
include("../profiling/profile.jl")
