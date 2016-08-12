######################################################
#
# Using adjacent frequencies for noise estimate
#
######################################################


"""
    beamformer_lcmv(s::SSR, l::Leadfield)

Linearly constrained minimum variance (LCMV) beamformer for epoched data.

NAI is ratio between frequeny of interest (`foi`) and specified control (`foi + noise_delta`) frequency.


### Input

* s = stimulus condition SSR data with epochs pre calculated
* l = leadfield information

* foi = frequency of interest for cross power spectral density calculations
* fs = sample rate
* n_epochs = number of epochs to average down to with aim of reducing noise

"""
function beamformer_lcmv(s::SSR, l::Leadfield;
                         foi::Real=modulationrate(s), fs::Real=samplingrate(s), n_epochs::Int=0,
                         freq_pm::Real = 0.5, noise_delta::Real=4.0, bilateral::Real = 15, kwargs...)

    Logging.info("Performing LCMV beamforming on signal data with $(foi + noise_delta) Hz as reference")

    if !haskey(s.processing, "epochs")
          Logging.critical("Epochs not calculated")
    end

    if n_epochs > 0
        s.processing["epochs"] = reduce_epochs(s.processing["epochs"], n_epochs)
    end

    l = match_leadfield(l, s)

    C = cross_spectral_density(s.processing["epochs"], foi - freq_pm, foi + freq_pm, fs; kwargs...)
    Q = cross_spectral_density(s.processing["epochs"], foi + noise_delta - freq_pm, foi + noise_delta + freq_pm, fs, ignore = foi; kwargs...)

    @assert size(C) == size(Q)
    @assert C != Q

    Logging.debug("Covariance matrices calculated and of size $(size(Q))")

    V, N, NAI = beamformer_lcmv(C, Q, l.L, l.x, l.y, l.z, bilateral; kwargs...)

    VolumeImage(vec(NAI), "NAI", l.x, l.y, l.z, ones(size(vec(NAI))), "LCMV", Dict(), "Talairach")
end
