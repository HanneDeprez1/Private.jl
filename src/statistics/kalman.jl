using TimeModels
using EEG

function kalman_filt(s::SSR; freq_of_interest::Union(Real, AbstractArray) = modulationrate(s), ID::String = "",
            model_type::Function=acoustic_model, proc_noise_cov=1e-10, obs_noise_cov=cov(s.data), state_noise_cov=0.001, x0=[0.0, 0.0],
            results_key::String="statistics", amp_est_start::Int = int(round(samplingrate(s) * 0.10)), kwargs...)

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

            amp = median(model_amplitude( kalman_filter( s.data[:, elec_idx], model) )[amp_est_start:end])

            result = DataFrame( ID                  = vec([ID]),
                                Channel             = vec([elec]),
                                ModulationRate      = modulationrate(s),
                                AnalysisType        = vec(["Kalman"]),
                                AnalysisFrequency   = freq,
                                SignalAmplitude     = amp)

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


function model_amplitude{T <: Number}(filte::KalmanFiltered{T})

    sqrt(filte.filtered[:, 1].^2 .+ filte.filtered[:, 2].^2)
end




"""
Build a basic SSR model with single sinusoid at the modulation rate.

Returns a state space model.

"""
function acoustic_model{T}(modulation_rate::T, sample_rate::T, v::T, w::T, p::T, x0::Vector{T}; kwargs...)

    Δt = 1 / sample_rate
    ω = 2 * pi * modulation_rate
    n = length(x0)

    ## Process transition and noise covariance

    F = eye(n)
    V = v * eye(n)  # We expect little variation in the state matrix

    ## Observation and noise covariance

    function G1(k);  cos(ω * Δt * k); end
    function G2(k); -sin(ω * Δt * k); end
    G = reshape([G1, G2], 1, n)
    W = diagm([w])   # Noise Covariance

    ## Inital guesses at state and error covariance
    # x0 set in function
    P0 = p * eye(n)

    StateSpaceModel(F, V, G, W, x0, P0)
end


"""
Build a SSR model with single sinusoid at the modulation rate and artifacts at the carrier rate.

Returns a state space model.

"""
function artifact_model{T}(modulation_rate::T, sample_rate::T, v::T, w::T, p::T, x0::Vector{T}; kwargs...)

    # Put in check that fa was passed in and v2

    Δt = 1 / sample_rate
    ω = 2 * pi * modulation_rate
    n = length(x0)

    ## Process transition and noise covariance

    F = eye(n)
    V = v * eye(n)  # We expect little variation in the state matrix

    ## Observation and noise covariance

    function G1(k);  cos(ω * Δt * k); end
    function G2(k); -sin(ω * Δt * k); end
    function G3(k);  triangle_artifact(Δt * k, fa; kwargs...); end
    G = reshape([G1, G2], 1, n)
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



"""
Polynomial smoothing with the Savitsky Golay filters

 https://github.com/blakejohnson/Qlab.jl/blob/master/src/SavitskyGolay.jl

 Sources
 ---------
 Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf
 Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay
"""
function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

    #Some error checking
    @assert isodd(windowSize) "Window size must be an odd integer."
    @assert polyOrder < windowSize "Polynomial order must me less than window size."

    halfWindow = int((windowSize-1)/2)

    #Setup the S matrix of basis vectors. 
    S = zeros(windowSize, polyOrder+1)
    for ct = 0:polyOrder
        S[:,ct+1] = [-halfWindow:halfWindow].^(ct)
    end

    #Compute the filter coefficients for all orders
    #From the scipy code it seems pinv(S) and taking rows should be enough
    G = S*pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv);

    #Pad the signal with the endpoints and convolve with filter
    paddedX = [x[1]*ones(halfWindow), x, x[end]*ones(halfWindow)]
    y = conv(filterCoeffs[end:-1:1], paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]

end

using LsqFit
using Loess

"""
Find trend for SSR type
"""
function trend(a::SSR; secs=5, order =3)

    ls = size(a.data, 1)

    windowSize = int(round(8192 * secs))
    if iseven(windowSize)
        windowSize += 1
    end

    savitsky_golay(vec(a.data), windowSize, order)

end

"""
Remove trend from SSR type
"""
function detrend(a::SSR; kwargs...)

    a.data = a.data - trend(a; kwargs...)

    return a
end
