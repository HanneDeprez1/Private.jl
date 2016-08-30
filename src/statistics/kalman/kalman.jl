using TimeModels, DSP

"""
    kalman_filter(s::SSR; kwargs...)

Use a kalman filter to estimate signal amplitude

#### Arguments

* 'freq_of_interest': Frequencies to analyse (modulation rate)
* 'model_type': Type of model to use (single channel ASSR)
* 'ID': Value to store as ID (empty "")
* 'results_key': Where to store results in processing dictionary ("statistics")
* 'proc_noise_cov': Process noise covariance (1e-10)
* 'obs_noise_cov': Observation noise covariance (covariance of signal)
* 'error_cov': Error covariance, measure of accuracy of state estimate (1.0)
* 'x0': Initial state estimates

"""
function TimeModels.kalman_filter(s::SSR;
            freq_of_interest::Union{Real, AbstractArray} = modulationrate(s),
            model_type::Function = acoustic_model, reduction_method::Function = mean,
            ID::AbstractString = "", results_key::AbstractString = "statistics",
            proc_noise_cov = 1e-10, obs_noise_cov = cov(s.data[find(~isnan(s.data[:, 1])), :]),
            error_cov = 1.0, x0=[0.0, 0.0], kwargs...)

    Logging.info("Running single channel Kalman filter on SSR data for $(size)")

    # Convert data to 64 bit type
    s.data = convert(Array{Float64, 2}, s.data)
    obs_noise_cov = convert(Array{typeof(proc_noise_cov)}, obs_noise_cov)

    for freq in freq_of_interest
        for elec in channelnames(s)

            elec_idx = findfirst(channelnames(s), elec)

            debug("Processing channel $elec in index $elec_idx")

            model = model_type(freq, samplingrate(s), proc_noise_cov, obs_noise_cov[elec_idx, elec_idx], error_cov, x0)

            f = kalman_smooth(s.data[:, elec_idx], model)

            amp = reduction_method(model_amplitude(f, s))
            noi = std(model_amplitude_time(f, s))
            pha = reduction_method(DSP.unwrap(model_phase(f), range = pi))

            result = DataFrame( ID                  = vec([ID]),
                                Channel             = vec([elec]),
                                ModulationRate      = modulationrate(s),
                                AnalysisType        = vec(["Kalman"]),
                                AnalysisFrequency   = freq,
                                SignalAmplitude     = amp,
                                SignalPhase         = pha,
                                NoiseAmplitude      = noi,
                                SNRdB               = 10 * log10((amp^2) / (noi^2)),
                                StateVariable1      = f.smoothed[end, 1],
                                StateVariable2      = f.smoothed[end, 2],
                                ProcessCov          = proc_noise_cov,
                                ObservationCov      = obs_noise_cov[elec_idx, elec_idx],
                                ErrorCov            = error_cov)

            result = add_dataframe_static_rows(result, kwargs)

            if haskey(s.processing, results_key)
                s.processing[results_key] = vcat(s.processing[results_key], result)
            else
                s.processing[results_key] = result
            end
        end
    end

    return s
end


"""
Return the amplitude for the SSR state space model
"""
function model_amplitude{T <: Number}(filte::Union{KalmanSmoothed{T}, TimeModels.KalmanSmoothedMinimal{T}}, s::SSR)

    amp = filte.smoothed[find(~isnan(s.data[:, 1])), :]
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude{T <: Number}(filte::Union{KalmanSmoothed{T}, TimeModels.KalmanSmoothedMinimal{T}})

    amp = filte.smoothed
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude_time{T <: Number}(filte::Union{KalmanSmoothed{T}, TimeModels.KalmanSmoothedMinimal{T}}, s::SSR)

    amp = filte.smoothed[find(~isnan(s.data[:, 1])), :]
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude_time{T <: Number}(filte::Union{KalmanSmoothed{T}, TimeModels.KalmanSmoothedMinimal{T}})

    amp = filte.smoothed
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_phase{T <: Number}(filte::Union{KalmanSmoothed{T}, TimeModels.KalmanSmoothedMinimal{T}})

    pha = filte.smoothed
    pha = atan(pha[:, 2] ./ pha[:, 1])
end

function model_phase{T <: Number}(filte::KalmanFiltered{T})

    pha = filte.filtered
    pha = atan(pha[:, 2] ./ pha[:, 1])
end


"""
Build a basic SSR model with single sinusoid at the modulation rate.

Returns a state space model.

#### Parameters

* 'modulation_rate': Modulation rate of sinusoid to match
* 'sample_rate': Sample rate of signal
* 'v': Process noise covariance (1e-10)
* 'w': Observation noise covariance (covariance of signal)
* 'p': Error covariance, measure of accuracy of state estimate (1.0)
* 'x0': Initial state estimates
* 'num_sensors': Number of sensors

"""
function acoustic_model{T}(modulation_rate::T, sample_rate::T, v::T, w::T, p::T, x0::Vector{T};
            num_sensors::Int = 1, kwargs...)

    num_states = length(x0)                 # Number of state variables
    V = v * eye(num_states)                 # Process noise covariance
    W = w * eye(num_sensors)                # Observation noise covariance
    P0 = p * eye(num_states)                # Initial noice covariance

    acoustic_model(num_sensors, modulation_rate, sample_rate, V, W, P0, x0; kwargs...)
end

function acoustic_model{T}(num_sensors::Int, modulation_rate::T, sample_rate::T,
            V::Array{T, 2}, W::Array{T, 2}, P0::Array{T, 2}, x0::Array{T, 1}; kwargs...)

    Δt = 1 / sample_rate
    ω = 2 * pi * modulation_rate
    num_states = length(x0)

    ## Process transition
    F = eye(num_states)

    # G = Array(Function, num_sensors, num_states)
    #
    # ## Observation
    # for sensor in 1:num_sensors
    #     G[sensor, 1] = k ->  cos(ω * Δt * k)
    #     G[sensor, 2] = k -> -sin(ω * Δt * k)
    # end

    function G(k::Int)
        [cos(ω * Δt * k); -sin(ω * Δt * k)]'
    end

    StateSpaceModel(_->F, _->V, G, _->W, x0, P0)
end


"""
Build a SSR model with single sinusoid at the modulation rate and artifacts at the carrier rate.

Returns a state space model.

"""
function artifact_model{T}(modulation_rate::T, sample_rate::T, v::Vector{T}, w::T, p::Vector{T}, x0::Vector{T}, fa::T; kwargs...)

    # Put in check that fa was passed in and v2

    Δt = 1 / sample_rate
    ω = 2 * pi * modulation_rate
    n = length(x0)

    ## Process transition and noise covariance

    F = eye(n)
    V = diagm(v)  # We expect little variation in the state matrix

    ## Observation and noise covariance

    function G1(k);  cos(ω * Δt * k); end
    function G2(k); -sin(ω * Δt * k); end
    function G3(k);  triangle_artifact(Δt * k, fa; kwargs...); end
    G = reshape([G1, G2, G3], 1, n)
    W = diagm([w])   # Noise Covariance

    ## Inital guesses at state and error covariance
    # x0 set in function
    P0 = diagm(p)

    StateSpaceModel(F, V, G, W, x0, P0)
end

function triangle_artifact(t::Float64, fa::Number; adj = 0.0, offset = 0.0, tol = 0.0001)

    dis1 = mod(-(t + offset), 1 / (fa+adj))
    dis2 = mod( (t + offset), 1 / (fa+adj))
    dis = min(dis1, dis2)

    if dis < tol
        x = 1.0 - abs(dis / tol)
    else
        x = 0.0
    end

    return x
end

function triangle_artifact(t::Array{Float64,1}, fa::Number; kwargs...)

    a = Array(Float64, size(t))

    for i in 1:length(t)
        a[i] = triangle_artifact(t[i], fa; kwargs...)
    end
    return a
end

"""
    set_nans(a)

Set samples to nan that are likely to be artifacts.

Set values to be excluded from analysis to NaN.
Usually in ASSR analysis you would remove the epochs, but for Kalman filter the timing would then be lost.
Instead we set the value to NaN and the filter does not updates its state estimates over these values, but does propogate the internal system.
The artifact ± `boundary` are set to NaN. The first and last 5 samples are not altered.

Currently only works on a single channel.
"""
function set_nans(a::SSR; boundary::Int = round(Int, 8192/2), perc=0.1)

    @assert size(a.data, 2) == 1 "Set NaNs only works for single channel SSR"

    info("Removing artifacts (NaN method)")
    debug("Input data consists of $(round(sum(isnan(a.data)) / length(a.data), 3))% NaNs")

    # We want to maintain the first and last 5 samples for Kalman filtering
    # Save them here and replace them at the end in case they were removed
    start_data = a.data[1:5, :]
    end_data = a.data[end-5:end, :]

    while (sum(isnan(a.data)) / length(a.data)) <= perc

        idx = findfirst(a.data, maximum(a.data))

        debug("Removing a segment of length $(2 * boundary) samples centered at sample $idx")

        new_min = maximum([1, idx - boundary])
        new_max = minimum([length(a.data), idx + boundary])

        a.data[new_min:new_max] = NaN
    end

    a.data[1:5, :] = start_data
    a.data[end-5:end, :] = end_data

    debug("Output data consists of $(round(sum(isnan(a.data)) / length(a.data), 3))% NaNs")

    return a
end
