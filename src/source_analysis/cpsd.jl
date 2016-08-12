"""
    cross_spectral_density{T <: AbstractFloat}(epochs::Array{T, 3}, fmin::Real, fmax::Real, fs::Real; fsum::Bool=false, ignore::Real=Inf, ignore_pm::Real=1.0, kwargs...)

Compute complex cross spectral density of epoch data

Cross spectral density averaged over frequencies of interest

### Input

* epochs: data shaped as epochs (samples x trials x channels)
* fmin: minimum frequency of interest
* max: maximum frequency of interest

### Implementation

Currently uses MNE python library.
Will change to Synchrony.jl when its stabilised.

"""
function cross_spectral_density{T <: AbstractFloat}(epochs::Array{T, 3}, fmin::Real, fmax::Real, fs::Real;
        fsum::Bool=false, ignore::Real=Inf, ignore_pm::Real=1.0, kwargs...)

    @pyimport mne as mne
    @pyimport mne.time_frequency as tf

    Logging.debug("Cross power spectral density between $fmin - $fmax Hz for signal of $fs Hz")

    # Convert from EEG.jl to MNE format for epochs
    epochs = permutedims(epochs, [2, 3, 1])                      # Change to trials x channels x samples
    names = AbstractString[string(i) for i = 1:size(epochs, 2)]  # Hack to put in fake names
    events = ones(Int, size(epochs, 1), 3)                       # Make all events the same to use everything
    events[:, 1] = 0:size(events, 1) - 1

    # Run MNE processing
    i = pycall(mne.create_info, PyObject, ch_names = vec(names), ch_types = vec(repmat(["eeg"], size(epochs, 2))), sfreq = fs)
    epochs = mne.EpochsArray(epochs, i, events, 0.0, verbose = false)

    csd = tf.compute_epochs_csd(epochs, fmin = fmin, fmax = fmax, verbose = false, fsum = fsum)

    keep_idx = find(abs(AbstractFloat[c[:frequencies][1] for c in csd] - ignore) .> ignore_pm)
    csd = csd[keep_idx];
    Logging.debug("Averaging CSD over frequencies $(AbstractFloat[c[:frequencies][1] for c in csd])")
    a = zeros(csd[1][:data])
    for i in 1:length(csd)
        a = a .+ csd[i][:data]
    end
    a = a ./ length(csd)

    return a
end

function cross_spectral_density(s::SSR; fmin::Real=modulationrate(s)-0.5, fmax::Real=modulationrate(s)+0.5,
        fs::Real=samplingrate(s), kwargs...)
    cross_spectral_density(s.processing["epochs"], fmin, fmax, fs; kwargs...)
end
