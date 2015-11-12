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

export
    import_leadfield,
    readLFT,
    readSRF,
    readLOC,
    readELP
include("types/Leadfield/Leadfield.jl")
include("types/Leadfield/import.jl")

end # module
