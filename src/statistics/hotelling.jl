using Distributions
using EEG
using DataFrames
using DSP
using Logging
using Docile
using Compat
using Gadfly


@doc """
Hotelling test on SSR data
Saves results in a.processing["hotelling"]
""" ->
function hotelling(s::SSR; freq_of_interest::Union{Real, AbstractArray} = modulationrate(s), ID::AbstractString = "",
    spectrum_type::Function = _hotelling_spectrum, data_type::AbstractString="epochs",
    fs::Number=samplingrate(s), results_key::AbstractString="statistics", kwargs...)

    # Calculate spectrum of each epoch
    # Do calculation here once, instead of in each low level call
    spectrum, frequencies = spectrum_type(s.processing[data_type], fs, freq_of_interest)
    spectrum  = compensate_for_filter(s.processing, spectrum, fs)

    for freq in freq_of_interest

        snrDb, phase, signal, noise, statistic = hotelling(spectrum, frequencies, freq)

        actual_freq = frequencies[_find_closest_number_idx(frequencies, freq)]

        result = DataFrame( ID                  = vec(repmat([ID], length(channelnames(s)), 1)),
                            Channel             = copy(channelnames(s)),
                            ModulationRate      = copy(modulationrate(s)),
                            AnalysisType        = "Hotelling",
                            AnalysisFrequency   = actual_freq,
                            SignalAmplitude     = vec(sqrt(signal)),
                            SignalPhase         = vec(phase),
                            NoiseAmplitude      = vec(sqrt(noise)),
                            SNRdB               = vec(snrDb),
                            Statistic           = vec(statistic)  )

        result = add_dataframe_static_rows(result, kwargs)

        if haskey(s.processing, results_key)
            s.processing[results_key] = vcat(s.processing[results_key], result)
        else
            s.processing[results_key] = result
        end

    end

    return s
end

# Backward compatibility
function hotelling(a::SSR, freq_of_interest::Union{Real, AbstractArray}; kwargs...)

    hotelling(a; freq_of_interest=freq_of_interest, kwargs...)
end


function hotelling{T <: AbstractFloat}(spectrum::Array{Complex{T},3}, frequencies::AbstractArray, freq_of_interest::Real)

    Logging.info("Calculating hotelling statistic on $(size(spectrum)[end]) channels at $freq_of_interest Hz with $(size(spectrum)[2]) epochs")

    idx  = _find_closest_number_idx(frequencies, freq_of_interest)
    for_stats = [real(spectrum[idx, :, :]) ; imag(spectrum[idx, :, :])]

    signal_amplitude = Array(T, size(spectrum, 3))
    signal_power     = Array(T, size(spectrum, 3))
    signal_phase     = Array(T, size(spectrum, 3))
    noise_power      = Array(T, size(spectrum, 3))
    snrDb            = Array(T, size(spectrum, 3))
    statistic        = Array(T, size(spectrum, 3))

    for c in 1:size(spectrum,3)

        # Signal metrics
        signal_amplitude[c] = abs(mean(spectrum[idx, :, c], 2))[1]
        signal_power[c]     = signal_amplitude[c]^2
        signal_phase[c]     = angle(mean(spectrum[idx, :, c], 2))[1]

        # Calculate noise
        n = [squeeze(spectrum[idx, :, c], 1); mean(squeeze(spectrum[idx, :, c], 1))]
        n = std(n)                      # Noise amplitude per epoch
        n = n / sqrt(size(spectrum, 2)) # Noise amplitude of recording

        # Save metrics
        noise_power[c] = real(n .^2)    # Noise power of recording
        snrDb[c]       = 10 * log10( signal_power[c] / noise_power[c] )
        statistic[c]   = _hotelling_T2_1sample(for_stats[:, :, c]')
    end

    return snrDb, signal_phase, signal_power, noise_power, statistic
end


function _hotelling_T2_1sample(data; corr::Number=1)

    n, p = size(data)

    m = mean(data, 1)

    if sum(cov(data)) == 0

        statistic = 1.0

    else

        T2 = (corr * n) * (m / cov(data)) * m'
        continuous_distribution = Chisq(p)
        statistic = ccdf(continuous_distribution, T2)[1]

    end

    return statistic
end



function _apes_spectrum{T <: AbstractFloat}(sweep::Array{T,3}, fs::Number, freq_of_interest)

    spectrum = apes(sweep, freq_of_interest / (fs / 2))

    return spectrum, freq_of_interest
end


# Calculate spectrum and return the associated frequencies
# Ignore the freq_of_interest variable as this is used for other spectrum estimation types
function _hotelling_spectrum{T <: AbstractFloat}(sweep::Array{T,3}, fs::Number, freq_of_interest)

    sweepLen = size(sweep)[1]
    
    p = plan_rfft(sweep, 1, flags = FFTW.MEASURE)

    spectrum = (2 / sweepLen) * (p * sweep)

    frequencies = linspace(0, 1, size(spectrum, 1)) * fs / 2

    return spectrum, frequencies
end


# For backwards compatibility. TODO remove
#function _hotelling_spectrum{T <: AbstractFloat}(sweep::Array{T,3})
#
#    sweepLen = size(sweep)[1]
#
#    (2 / sweepLen) * fft(sweep, 1)[1:sweepLen / 2 + 1, :, :]
#end


function plot_hotelling{T <: AbstractFloat}(spectrum::Array{Complex{T},3}, frequencies::AbstractArray, freq_of_interest::Real; c::Int=1, fig_name="hot.pdf")

    idx  = _find_closest_number_idx(frequencies, freq_of_interest)
    for_plots = [vec(real(spectrum[idx, :, :])) ; vec(imag(spectrum[idx, :, :]))]

    # Signal metrics
    signal_amplitude = abs(mean(spectrum[idx, :, c], 2))[1]
    signal_power     = signal_amplitude[c]^2
    signal_phase     = angle(mean(spectrum[idx, :, c], 2))[1]

    # Calculate noise
    n = [squeeze(spectrum[idx, :, c], 1), mean(squeeze(spectrum[idx, :, c], 1))]
    n = std(n)                      # Noise amplitude per epoch
    n = n / sqrt(size(spectrum, 2)) # Noise amplitude of recording

    # Save metrics
    noise_power = real(n .^2)    # Noise power of recording
    snrDb       = 10 * log10( signal_power / noise_power )

    d = layer(x=for_plots, Geom.histogram)
    s = layer(xintercept=[signal_amplitude], Geom.vline(color="black"))
    n = layer(xintercept=[-n, n] .+ signal_amplitude, Geom.vline(color="red"))
    p = plot(n, s, d,
        Guide.xlabel("Real and Imag FFT Bin (uV)"),
        Guide.title("SNR = $(round(snrDb[1], 3)) (dB)"))

    draw(PDF(fig_name, 26cm, 17cm), p)

    return p
end
