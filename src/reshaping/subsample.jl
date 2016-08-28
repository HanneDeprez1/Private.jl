function blank(a::SSR, blank_delay::AbstractFloat; valid_triggers::Int=-4,
    temptrigger_rate::AbstractFloat=1/a.processing["Carrier_Frequency"], temptrigger_rate_code::Int=22,
    temptrigger::Int=33, kwargs...)

    Logging.info("Blanking SSR after $blank_delay (s)")

    # Place trigger to indicate where electrical stimulation occurs (22)
    new_triggers = extra_triggers(a.triggers, valid_triggers, temptrigger_rate_code, temptrigger_rate, samplingrate(a))

    # Place trigger where to blank from
    new_triggers = extra_triggers(new_triggers, [valid_triggers; temptrigger_rate_code], temptrigger, blank_delay, samplingrate(a), max_inserted=1)

    # Create a time vector
    a.processing["time"] = collect(1 : size(a.data, 1)) / samplingrate(a)

    # Find points to blank from
    valid_trip_idx = find(new_triggers["Code"]-252 .== temptrigger)

    # Initial values
    previous_idx = 1
    previous_value = a.data[previous_idx, :]

    # For each point to sample
    for i in 1:length(valid_trip_idx)-1

        current_valid_idx = new_triggers["Index"][valid_trip_idx[i]]
        current_value = a.data[current_valid_idx, :]
        idxs = collect(previous_idx : 1 : current_valid_idx)

        for c = 1:size(a.data, 2) # For each channel

            a.data[idxs, c] = linspace(previous_value[c], current_value[c], length(idxs))
        end

        previous_value = current_value
        previous_idx = current_valid_idx

    end

    return a
end



function subsample(a::SSR; valid_triggers::Int=-4,
                   temptrigger_rate::AbstractFloat=1/a.processing["Carrier_Frequency"], temptrigger_rate_code::Int=22,
                   subsample_start_delay::Number=0.001, temptrigger_start_idx=33,
                   subsample_stop_delay::Number=0.0015, temptrigger_stop_idx=34,
                   plot::Bool=false, kwargs...)

    # Add a trigger (22) where each CI artifact starts
    # Add a trigger (33) where to start sampling the valid response
    # Add a trigger (34) where to stop sampling the valid response
    # Go through each (33) trigger and take mean between that and (34)
    # Place that mean value at the mean time point between them and save subsampled data to SSR type

    Logging.info("Subsampling SSR")
    Logging.debug("Pulse every $(round(temptrigger_rate, 8)) s")
    Logging.debug("Start sampling after $(round(subsample_start_delay, 8)) s")
    Logging.debug("Stop sampling after $(round(subsample_stop_delay, 8)) s")

    a.processing["time"] = collect(1:size(a.data,1)) / samplingrate(a)

    if plot
        warn("Plotting functionality has been depreciated")
    end

    # Place trigger to indicate where electrical stimulation occurs (22)
    new_triggers = extra_triggers(a.triggers, valid_triggers, temptrigger_rate_code, temptrigger_rate, samplingrate(a))

    # Place trigger to indicate where clean data starts (33)
    new_triggers = extra_triggers(new_triggers, [valid_triggers, temptrigger_rate_code], temptrigger_start_idx, subsample_start_delay, samplingrate(a), max_inserted=1)

    # Place trigger to indicate where clean data ends (34)
    new_triggers = extra_triggers(new_triggers, [valid_triggers, temptrigger_rate_code], temptrigger_stop_idx, subsample_stop_delay, samplingrate(a), max_inserted=1)

    # if length(unique(new_triggers["Index"])) != length(new_triggers["Index"])
    #     warn("Multiple triggers in the same place. Start and end of sample location might be too close together")
    # end

    # Run through all the 33s and interpolate between averaged valid values
    valid_trip_idx = find(new_triggers["Code"]-252 .== 33)

    previous_value = zeros(Float32, size(a.data,2)) # Not used, but will stop
    previous_idx = 1                                # lint errors

    cnt = 0
    for i in 1:length(valid_trip_idx)-1
        cnt += 1

        if cnt < 4
            Logging.info("Processing sample location $cnt with idx = $(valid_trip_idx[i]) and value $(new_triggers["Code"][valid_trip_idx[i]]-252)")
        end

        # Check the next index is a 34
        if new_triggers["Code"][valid_trip_idx[i]+1]-252 == temptrigger_stop_idx

            valid_range = new_triggers["Index"][valid_trip_idx[i]] : 1 : new_triggers["Index"][valid_trip_idx[i]+1]

            if i == 1
              Logging.debug("Averaging $(length(valid_range)) data points for subsample")
            end

            mean_value = mean(a.data[valid_range ,:], 1)

            mean_idx = int(round(mean(valid_range)))

            # Only interpolate from second valid point backwards
            if cnt > 1

                idxs = [previous_idx : 1 : mean_idx]

                for c = 1:size(a.data,2)

                    a.data[idxs, c] = linspace(previous_value[c], mean_value[c], length(idxs))

                end
            end

            previous_value = mean_value
            previous_idx = mean_idx

        else

            error("Expected next trigger $cnt to be end of valid range. But was $(new_triggers["Code"][valid_trip_idx[i]+1]-252)")
            error(new_triggers["Code"][1:10]-252)
            error(new_triggers["Index"][1:10])
        end
    end

    return a
end
