module Private

export
    hotelling,
    plot_hotelling
include("statistics/hotelling.jl")

export
    subsample
include("reshaping/subsample.jl")

end # module
