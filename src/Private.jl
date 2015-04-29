module Private

export
    hotelling,
    plot_hotelling
include("statistics/hotelling.jl")

export
    subsample,
    blank
include("reshaping/subsample.jl")

end # module
