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

    spectrum    = _hotelling_spectrum(a.processing["epochs"])
    frequencies = linspace(0, 1, int(size(spectrum, 1)))*float(a.sample_rate)/2

    # TODO: Account for multiple applied filters
    if haskey(a.processing, "filter1")
        used_filter         = a.processing["filter1"]
        filter_response     = freqz(used_filter, frequencies, float(a.sample_rate))
        for f = 1:length(filter_response)
            spectrum[f, :, :] = spectrum[f, :, :] ./ abs(filter_response[f])^2
        end
        debug("Accounted for filter response in hotelling test spectrum estimation")
    end

    for freq in freq_of_interest

        snrDb, phase, signal, noise, statistic = hotelling(spectrum, frequencies, freq, int(a.sample_rate))

        result = DataFrame(
                            ID                  = vec(repmat([ID], length(a.channel_names), 1)),
                            Channel             = copy(a.channel_names),
                            ModulationFrequency = copy(float(a.modulation_frequency)),
                            AnalysisType        = "hotelling",
                            AnalysisFrequency   = freq,
                            SignalPower         = vec(signal),
                            SignalPhase         = vec(phase),
                            NoisePower          = convert(Array{FloatingPoint}, noise),
                            SNRdB               = convert(Array{FloatingPoint}, snrDb),
                            Statistic           = vec(statistic)
                          )

        result = add_dataframe_static_rows(result, kwargs)

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

    idx         = _find_closest_number_idx(frequencies, freq_of_interest)

    bins = spectrum[idx, :, :]

    signal_amplitude = float(abs(squeeze(mean(bins,2), [1, 2])))
    signal_power     = float(signal_amplitude.^2)
    signal_phase     = float(angle(squeeze(mean(bins,2), [1, 2])))

    noise_power = FloatingPoint[]
    for i in 1:size(bins,3)
        d = [squeeze(bins[1, :, i], [1, 3]), mean(squeeze(bins[1, :, i], [1, 3]))]
        d = std(d)                      # epoch_noise_amplitude
        d = d / sqrt(size(spectrum, 2)) # recording_noise_amplitude
        d = real(d .^2)                 # recording_noise_power
        push!(noise_power, d)
    end

    snr   = signal_power ./ noise_power
    snrDb = 10*log10(snr)

    for_stats = [real(bins) , imag(bins)]
    statistic = [_hotelling_T2_1sample(for_stats[:, :, i]') for i = 1:size(for_stats, 3)]

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
