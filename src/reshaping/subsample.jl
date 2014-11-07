using Gadfly

function subsample(a::SSR; plot_channel::Int=1, subsample_start_delay::Number=0.001, subsample_stop_delay::Number=0.0015,
                   plot::Bool=false, plot_center::Number=2.9, plot_name::String="$(a.file_name)-subsample")

    # Add a trigger (22) where each CI artifact starts
    # Add a trigger (33) where to start sampling the valid response
    # Add a trigger (34) where to stop sampling the valid response
    # Go through each (33) trigger and take mean between that and (34)
    # Place that mean value at the mean time point between them and save subsampled data to SSR type

    info("Subsampling SSR")

    a.processing["time"] = [1:size(a.data,1)]/float(a.sample_rate)

    if plot

        l0 = layer( x=a.processing["time"], y=a.data[:,plot_channel], Geom.line)
        l1 = layer( xintercept=a.triggers["Index"]/float(a.sample_rate), Geom.vline(color="black"))

        p1 = Gadfly.plot(l0, l1,
            Scale.x_continuous(maxvalue=maximum(a.processing["time"])),
            #=Scale.x_continuous(minvalue=plot_center-6, maxvalue=plot_center+6),=#
            Guide.ylabel("Amplitude (uV)"),
            Guide.xlabel("Time (s)"),
            Guide.title("EEG $(a.channel_names[plot_channel]) Epochs")
            )
    end


    new_triggers = extra_triggers(a.triggers, 1, 22, 1/a.processing["Carrier_Frequency"], float(a.sample_rate))

    new_triggers = extra_triggers(new_triggers, [1, 22], 33, subsample_start_delay, float(a.sample_rate), max_inserted=1)

    new_triggers = extra_triggers(new_triggers, [1, 22], 34, subsample_stop_delay,  float(a.sample_rate), max_inserted=1)

    # Run through all the 33s and interpolate between averaged valid values
    valid_trip_idx = find(new_triggers["Code"]-252 .== 33)

    count = 0
    previous_value = 0     # Not used, but will stop
    previous_idx   = 1     # lint error
    for i in 1:length(valid_trip_idx)-1
        count += 1

        # Check the next index is a 34
        if new_triggers["Code"][valid_trip_idx[i]+1]-252 == 34

            valid_range = new_triggers["Index"][valid_trip_idx[i]] : 1 : new_triggers["Index"][valid_trip_idx[i]+1]

            mean_value  = mean(a.data[valid_range ,:], 1)

            mean_idx    = int(round(mean(valid_range)))

            # Only interpolate from second valid point backwards
            if count > 1

                idxs = [previous_idx : 1 : mean_idx]

                for c = 1:size(a.data,2)

                    a.data[idxs, c] = linspace(previous_value[c], mean_value[c], length(idxs))

                end
            end

            previous_value = mean_value
            previous_idx   = mean_idx

        else

            error("Expected next trigger $count to be end of valid range. But was $(new_triggers["Code"][valid_trip_idx[i]+1]-252)")
        end
    end


    if plot

        l2 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 22+252]/float(a.sample_rate),
            Geom.vline(color="red")
            )

        l3 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 33+252]/float(a.sample_rate),
            Geom.vline(color="green")
            )

        l4 = layer(
            xintercept=new_triggers["Index"][new_triggers["Code"].== 34+252]/float(a.sample_rate),
            Geom.vline(color="orange")
            )

        l5 = layer(
            x=a.processing["time"],
            y=a.data[:, plot_channel],
            Theme(default_color=color("purple")),
            Geom.line
            )

        l7 = layer(
            xintercept=a.triggers["Index"]/float(a.sample_rate),
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


