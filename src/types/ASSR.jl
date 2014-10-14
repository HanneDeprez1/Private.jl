
function hotelling(a::ASSR, freq_of_interest::Number; ID::String="", kwargs...)

    snrDb, phase, signal, noise, statistic =
        hotelling(a.processing["epochs"], float(a.modulation_frequency), int(a.sample_rate))

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
function hotelling(a::ASSR, freq_of_interest::Array; kwargs...)

    for f = freq_of_interest; a = hotelling(a, f; kwargs...); end; return a
end


# If no frequency of interest is specified then use the modulation frequency
function hotelling(a::ASSR; kwargs...)

    hotelling(a, float(a.modulation_frequency); kwargs...)
end
