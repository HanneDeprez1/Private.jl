using Distributions
using EEG
using DataFrames
using DSP
using Logging
using Docile
using Compat
@docstrings


@doc """
Hotelling test on SSR data
Saves results in a.processing["hotelling#"]
""" ->
function hotelling(a::SSR; freq_of_interest::Union(Real, AbstractArray)=float(a.modulation_frequency), ID::String="", kwargs...)

    # Calculate spectrum of each epoch
    # Do calculation here once, instead of in each low level call
    spectrum    = _hotelling_spectrum(a.processing["epochs"])
    spectrum    = compensate_for_filter(a.processing, spectrum, float(a.sample_rate))
    frequencies = linspace(0, 1, int(size(spectrum, 1)))*float(a.sample_rate)/2

    for freq in freq_of_interest

        snrDb, phase, signal, noise, statistic = hotelling(spectrum, frequencies, freq, float(a.sample_rate))

        result = DataFrame( ID                  = vec(repmat([ID], length(a.channel_names), 1)),
                            Channel             = copy(a.channel_names),
                            ModulationFrequency = copy(float(a.modulation_frequency)),
                            AnalysisType        = "hotelling",
                            AnalysisFrequency   = freq,
                            SignalPower         = vec(signal),
                            SignalPhase         = vec(phase),
                            NoisePower          = vec(noise),
                            SNRdB               = vec(snrDb),
                            Statistic           = vec(statistic)  )

        result   = add_dataframe_static_rows(result, kwargs)
        key_name = new_processing_key(a.processing, "hotelling")
        merge!(a.processing, @compat Dict(key_name => result) )
    end

    return a
end

# Backward compatibility
function hotelling(a::SSR, freq_of_interest::Union(Real, AbstractArray); kwargs...)

    hotelling(a; freq_of_interest=freq_of_interest, kwargs...)
end


function hotelling(epochs::Union(Array{Float64, 3}, Array{Float32, 3}), args...; kwargs...)

    spectrum    = _hotelling_spectrum(a.processing["epochs"])
    frequencies = linspace(0, 1, int(size(spectrum, 1)))*float(a.sample_rate)/2

    hotelling(spectrum, frequencies, args...; kwargs...)
end

function hotelling(spectrum::Union(Array{Complex{Float64},3}, Array{Complex{Float32},3}, Array{Complex{FloatingPoint},3}),
                   frequencies::AbstractArray, freq_of_interest::Real, fs::Real)

    info("Calculating hotelling statistic on $(size(spectrum)[end]) channels at $freq_of_interest Hz with $(size(spectrum)[2]) epochs")

    idx  = _find_closest_number_idx(frequencies, freq_of_interest)
    for_stats = [real(spectrum[idx, :, :]) , imag(spectrum[idx, :, :])]

    signal_amplitude = Array(FloatingPoint, size(spectrum, 3))
    signal_power     = Array(FloatingPoint, size(spectrum, 3))
    signal_phase     = Array(FloatingPoint, size(spectrum, 3))
    noise_power      = Array(FloatingPoint, size(spectrum, 3))
    snrDb            = Array(FloatingPoint, size(spectrum, 3))
    statistic        = Array(FloatingPoint, size(spectrum, 3))

    for c in 1:size(spectrum,3)

        # Signal metrics
        signal_amplitude[c] = abs(mean(spectrum[idx, :, c], 2))[1]
        signal_power[c]     = signal_amplitude[c]^2
        signal_phase[c]     = angle(mean(spectrum[idx, :, c], 2))[1]

        # Calculate noise
        n = [squeeze(spectrum[idx, :, c], [1, 3]), mean(squeeze(spectrum[idx, :, c], [1, 3]))]
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


# Calculates the spectrum for hotelling
function _hotelling_spectrum(sweep::Array{Float64,3}; ref::Int=0)
    # First dimension is samples, second dimension if existing is channels

    sweepLen      = size(sweep)[1]

    # Calculate amplitude sweepe at each frequency along first dimension
    fftSweep    = 2 / sweepLen * fft(sweep, 1)

    spectrum    = fftSweep[1:sweepLen / 2 + 1, :, :]

    if ref > 0
        refspec = spectrum[:,ref]
        for i = 1:size(spectrum)[2]
            spectrum[:,i] = spectrum[:,i] - refspec
        end
    end

    return spectrum
end
