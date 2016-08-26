function blank(a::SSR, blank_delay::AbstractFloat; valid_triggers::Int=-4,
    temptrigger_rate::AbstractFloat=1/a.processing["Carrier_Frequency"], temptrigger_rate_code::Int=22,
    temptrigger::Int=33)

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
                   plot::Bool=false, plot_channel::Int=1, plot_center::Number=0.2, plot_name::AbstractString="$(a.file_name)-subsample")

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

        l0 = layer( x=a.processing["time"], y=a.data[:,plot_channel], Geom.line)
        l1 = layer( xintercept=a.triggers["Index"]/samplingrate(a), Geom.vline(color="black"))

        p1 = Gadfly.plot(l0, l1,
            Scale.x_continuous(maxvalue=maximum(a.processing["time"])),
            #=Scale.x_continuous(minvalue=plot_center-6, maxvalue=plot_center+6),=#
            Guide.ylabel("Amplitude (uV)"),
            Guide.xlabel("Time (s)"),
            Guide.title("EEG $(a.channel_names[plot_channel]) Epochs")
            )
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


    if plot

        l2 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 22+252]/samplingrate(a),
            Geom.vline(color="red")
            )

        l3 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 33+252]/samplingrate(a),
            Geom.vline(color="green")
            )

        l4 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 34+252]/samplingrate(a),
            Geom.vline(color="orange")
            )

        l5 = layer(
            x=a.processing["time"],
            y=a.data[:, plot_channel],
            Theme(default_color=color("purple")),
            Geom.line
            )

        l7 = layer(
            xintercept=a.triggers["Index"]/samplingrate(a),
            Geom.vline(color="black")
            )

        p2 = Gadfly.plot(
            l2,
            l1,
            l0,
            Scale.x_continuous(minvalue=plot_center-0.05, maxvalue=plot_center+0.05),
            Guide.ylabel("Amplitude (uV)"),
            Guide.xlabel("Time (s)"),
            Guide.title("EEG $(a.channel_names[plot_channel]) Pulses Per Second")
            )

        p3 = Gadfly.plot(
            l1,
            l2,
            l3,
            l4,
            l0,
            l5,
            Scale.x_continuous(minvalue=plot_center-0.01, maxvalue=plot_center+0.01),
            Guide.ylabel("Amplitude (uV)"),
            Guide.xlabel("Time (s)"),
            Guide.title("EEG $(a.channel_names[plot_channel]) Sampled Intervals")
            )

        p4 = Gadfly.plot(
            l7,
            l5,
            Scale.x_continuous(maxvalue=maximum(a.processing["time"])),
            #=Scale.x_continuous(minvalue=plot_center-6, maxvalue=plot_center+6),=#
            Guide.ylabel("Amplitude (uV)"),
            Guide.xlabel("Time (s)"),
            Guide.title("EEG $(a.channel_names[plot_channel]) Sub Sampled Data")
            )

        #=draw(SVGJS("$plot_name-3.js.svg", 14inch, 12inch), vstack(p1,p2,p3,p4))=#
        #=draw(PNG("$plot_name.png", 14inch, 12inch), vstack(p1,p2,p3,p4))=#
        draw(PNG("$plot_name.png", 11.02inch, 8.27inch), vstack(p1,p2,p3,p4))
        draw(PS("$plot_name.ps", 11.02inch, 8.27inch), vstack(p1,p2,p3,p4))

    end

    return a
end
