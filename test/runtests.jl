using Private
using Base.Test
using Logging
using Lint

Logging.configure(level=DEBUG)

include("source_analysis/beamformer.jl")
include("statistics/hotelling.jl")
include("statistics/kalman.jl")
include("reshaping/subsample.jl")
include("types/Leadfield/Leadfield.jl")

@printf("\n\n\n All tests passed. Now checking lint\n\n\n")

lintpkg( "Private" )
