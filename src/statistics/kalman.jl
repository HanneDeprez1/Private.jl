function TimeModels.kalman_filter(s::SSR; freq_of_interest::Union(Real, AbstractArray) = modulationrate(s),
            model_type::Function = acoustic_model, ID::String = "", x0=[0.0, 0.0], results_key::String = "statistics",
            proc_noise_cov = 1e-10, obs_noise_cov = cov(s.data), state_noise_cov = 0.001,
            amp_est_start::Int = int(round(samplingrate(s) * 0.10)), kwargs...)

    info("Running kalman filter analysis on SSR data")

    # Convert data to 64 bit type
    s.data = convert(Array{Float64, 2}, s.data)

    # Ensure parameters are the same type
    obs_noise_cov = convert(Array{typeof(proc_noise_cov)}, obs_noise_cov)

    for freq in freq_of_interest
        for elec in s.channel_names

            elec_idx = findfirst(s.channel_names, elec)

            debug("Processing channel $elec in index $elec_idx")

            model = model_type(freq, samplingrate(s), proc_noise_cov, obs_noise_cov[elec_idx, elec_idx], state_noise_cov, x0)

            f = kalman_smooth( s.data[:, elec_idx], model)

            amp = median(model_amplitude(f)[amp_est_start:end])
            noi = std(model_amplitude(f)[amp_est_start:end])

            result = DataFrame( ID                  = vec([ID]),
                                Channel             = vec([elec]),
                                ModulationRate      = modulationrate(s),
                                AnalysisType        = vec(["Kalman"]),
                                AnalysisFrequency   = freq,
                                SignalAmplitude     = amp,
                                NoiseAmplitude      = noi,
                                StateVariable1      = f.filtered[end, 1],
                                StateVariable2      = f.filtered[end, 2])

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
function model_amplitude{T <: Number}(filte::KalmanSmoothed{T}, s::SSR)

    amp = filte.smoothed[find(~isnan(s.data)), :]
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude{T <: Number}(filte::KalmanFiltered{T}, s::SSR)

    amp = filte.filtered[find(~isnan(s.data)), :]
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude{T <: Number}(filte::KalmanSmoothed{T})

    amp = filte.smoothed
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end

function model_amplitude{T <: Number}(filte::KalmanFiltered{T})

    amp = filte.filtered
    amp = sqrt(amp[:, 1].^2 .+ amp[:, 2].^2)
end


"""
Build a basic SSR model with single sinusoid at the modulation rate.

Returns a state space model.

"""
function acoustic_model{T}(modulation_rate::T, sample_rate::T, v::T, w::T, p::T, x0::Vector{T}; kwargs...)

    n = length(x0)  # Number of state variables
    V = v * eye(n)  # Process noise covariance
    W = diagm([w])  # Observation noise covariance
    P0 = p * eye(n) # Initial

    acoustic_model(modulation_rate, sample_rate, V, W, P0, x0; kwargs...)
end

function acoustic_model{T}(modulation_rate::T, sample_rate::T, V::Array{T, 2}, W::Array{T, 2}, P0::Array{T, 2}, x0::Array{T, 1}; kwargs...)

    Δt = 1 / sample_rate
    ω = 2 * pi * modulation_rate
    n = 2

    ## Process transition
    F = eye(2)

    ## Observation
    function G1(k);  cos(ω * Δt * k); end
    function G2(k); -sin(ω * Δt * k); end
    G = reshape([G1, G2], 1, n)

    StateSpaceModel(F, V, G, W, x0, P0)
end


"""
Build a SSR model with single sinusoid at the modulation rate and artifacts at the carrier rate.

Returns a state space model.

"""
function artifact_model{T}(modulation_rate::T, sample_rate::T, v::Vector{T}, w::T, p::T, x0::Vector{T}; kwargs...)

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
    G = reshape([G1, G2, G3, 1], 1, n)
    W = diagm([w])   # Noise Covariance

    ## Inital guesses at state and error covariance
    # x0 set in function
    P0 = p * eye(n)

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



