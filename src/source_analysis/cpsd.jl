###############################
#
# Cross power spectral density
#
###############################


function cross_power_spectral_density{T <: AbstractFloat}(epochs::Array{T, 3}, fmin::Number, fmax::Number, fs::Number)

    @pyimport mne as mne
    @pyimport mne.time_frequency as tf

    # Convert from EEG.jl to MNE format for epochs
    epochs = permutedims(epochs, [2, 3, 1])                      # Change to trials x channels x samples
    names = AbstractString[string(i) for i = 1:size(epochs, 2)]  # Hack to put in fake names
    events = ones(Int, size(epochs, 1), 3)                       # Make all events the same to use everything

    # Run MNE processing
    i = pycall(mne.create_info, PyObject, ch_names = vec(names), ch_types = vec(repmat(["eeg"], size(epochs, 2))), sfreq = fs)
    epochs = mne.EpochsArray(epochs, i, events, 0.0)
    csd = tf.compute_epochs_csd(epochs, fmin = fmin, fmax = fmax)
    csd = csd[:data]
    cpsd = abs(csd .* csd)
end

