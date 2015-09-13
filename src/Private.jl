module Private

using EEG         # General
using TimeModels  # Kalman
using LsqFit      # Detrend
using Loess       # Detrend

export
    kalman_filter,
    set_nans,
	model_amplitude
include("statistics/kalman.jl")

export
    hotelling,
    plot_hotelling
include("statistics/hotelling.jl")

export
    subsample,
    blank
include("reshaping/subsample.jl")

export
    detrend,
    trend
include("reshaping/detrend.jl")

end # module
