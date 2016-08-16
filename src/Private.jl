module Private

using EEG         # General
using TimeModels  # Kalman
using LsqFit      # Detrend
using Loess       # Detrend
using PyCall      # Cross Spectral Density
using SIUnits, SIUnits.ShortUnits
using ProgressMeter
using DataFrames

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
    readLFT,
    readSRF,
    readLOC,
    readELP
include("read_write/elp.jl")
include("read_write/srf.jl")
include("read_write/loc.jl")
include("read_write/lft.jl")

export
    import_leadfield,
    calculate_specific_leadfield,
    calculate_specific_leadfield2
include("types/Leadfield/Leadfield.jl")
include("types/Leadfield/processing.jl")

export
    beamformer_lcmv,
    cross_spectral_density,
    reduce_epochs,
    retain_svd,
    correct_midline
include("source_analysis/lcmv.jl")
include("source_analysis/cpsd.jl")
include("source_analysis/misc.jl")


end # module
