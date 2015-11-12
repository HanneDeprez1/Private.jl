using Private
using Base.Test
using Logging
using Lint

Logging.configure(level=DEBUG)

include("statistics/hotelling.jl")
include("statistics/kalman.jl")
include("reshaping/subsample.jl")
include("types/Leadfield/Leadfield.jl")
#= include("../profiling/profile.jl") =#
lintpkg( "Private" )
