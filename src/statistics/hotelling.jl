using Distributions
using EEG
using DataFrames
using DSP
using Logging
using Docile
@docstrings


@doc """
Hotelling test on SSR data
Saves results in a.processing["hotelling#"]
""" ->
function hotelling(a::SSR, freq_of_interest::Number; ID::String="", kwargs...)

    # TODO: Account for multiple applied filters
    if haskey(a.processing, "filter1")
        used_filter = a.processing["filter1"]
    else
        used_filter = nothing
    end

    snrDb, phase, signal, noise, statistic =
        hotelling(a.processing["epochs"], freq_of_interest, int(a.sample_rate), used_filter=used_filter; kwargs...)

    result = DataFrame(
                        ID                  = vec(repmat([ID], length(a.channel_names), 1)),
                        Channel             = copy(a.channel_names),
                        ModulationFrequency = copy(float(a.modulation_frequency)),
                        AnalysisType        = "hotelling",
                        AnalysisFrequency   = freq_of_interest,
                        SignalPower         = vec(signal),
                        SignalPhase         = vec(phase),
                        NoisePower          = convert(Array{FloatingPoint}, noise),
                        SNRdB               = convert(Array{FloatingPoint}, snrDb),
                        Statistic           = vec(statistic)
                      )

    result = add_dataframe_static_rows(result, kwargs)

    key_name = new_processing_key(a.processing, "hotelling")
    merge!(a.processing, [key_name => result])

    return a
end


# If more than one frequency of interest is specified then run for all
function hotelling(a::SSR, freq_of_interest::Array; kwargs...)

    for f = freq_of_interest; a = hotelling(a, f; kwargs...); end; return a
end


# If no frequency of interest is specified then use the modulation frequency
function hotelling(a::SSR; kwargs...)

    hotelling(a, float(a.modulation_frequency); kwargs...)
end




function hotelling(epochs::Array, freq_of_interest::Number, fs::Number;
                   used_filter=nothing, kwargs...)


    spectrum    = _hotelling_spectrum(epochs)
    frequencies = linspace(0, 1, int(size(epochs,1) / 2 + 1))*fs/2
    idx         = _find_closest_number_idx(frequencies, freq_of_interest)

    info("Calculating hotelling statistic on $(size(epochs)[end]) channels at $freq_of_interest Hz with $(size(epochs)[2]) epochs")

    # Compensate for filter response
    if !(used_filter == nothing)
        filter_response     = freqz(used_filter, frequencies, fs)
        filter_compensation = [abs(f)^2 for f = filter_response]
        spectrum            = spectrum ./ filter_compensation
        debug("Accounted for filter response in F test spectrum estimation")
    end

    bins = spectrum[idx, :, :]

    signal_amplitude = abs(squeeze(mean(bins,2), [1, 2]))
    signal_power     = signal_amplitude.^2
    signal_phase     = angle(squeeze(mean(bins,2), [1, 2]))


    noise_power = FloatingPoint[]
    for i in 1:size(bins,3)
        d = [squeeze(bins[1, :, i], [1, 3]), mean(squeeze(bins[1, :, i], [1, 3]))]
        d = std(d)                      # epoch_noise_amplitude
        d = d / sqrt(size(epochs,2))    # recording_noise_amplitude
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



# Calculates the spectrum for ftest and plotting
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


function second_moment(v)

    sqrt(sum((v - mean(v)).^2) / (length(v)))

end
