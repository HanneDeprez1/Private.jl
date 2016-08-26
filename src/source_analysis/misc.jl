###############################
#
# Retain SVD components
#
###############################

"""
Retain eigenvectors that represent up to `k` percent of power
"""
function retain_svd{T <: AbstractFloat}(A::Array{T, 2}, k::T=0.9)

    ss = svdfact(A)

    pw = cumsum(ss.S ./ sum(ss.S))

    keep = maximum(find(pw .< k)) + 1

    ss = ss.U[:, 1:keep]

    return ss, pw[keep]
end


###############################
#
# Average the epochs
#
###############################

"""
    reduce_epochs(Array{T, 3}, new_num_epochs::Int=30)

Average epochs down to specified number of epochs.
"""
function reduce_epochs{T <: AbstractFloat}(a::Array{T, 3}, new_num_epochs::Int=30)
    if new_num_epochs < size(a, 2)
        ep_per_av = floor(size(a, 2) / new_num_epochs)
        new = zeros(size(a, 1), new_num_epochs, size(a, 3))
        for i in 1:new_num_epochs-1
            new[:, i, :] = mean(a[:, 1+((i-1)*ep_per_av):ep_per_av+((i-1)*ep_per_av), :], 2)
        end
        new[:, new_num_epochs, :] = mean(a[:, 1+((new_num_epochs-1)*ep_per_av):end, :], 2)
        return new
    else
        return a
    end
end

function reduce_epochs(s::SSR,  new_num_epochs::Int=30)

    s.processing["epochs"] = reduce_epochs(s.processing["epochs"], new_num_epochs)
    return s
end

###############################
#
# Correct midline for bilateral
#
###############################


"""
Correct midline estimates from bilateral beamformer

Region supression for bilateral beamformers supresses itself along the midline. This causes innacurate estimates
around 0 on the x axis. To compensate for this take the average of adjacent valid locations.

#### Input

* vi: Volume image
* pm: Locations plus or minus the midline in mm to correct
"""
function correct_midline(v::VolumeImage; pm::Real=0.005, units=Meter)

    Logging.info("Correcting volume image midline errors caused by bilateral region supression")
    Logging.debug("Supressing region ± $(pm * units)")

    x = AbstractFloat[xi / (1 * units) for xi in v.x]
    midline_idxs = abs(x)
    midline_idxs = midline_idxs .<= pm

    valid_idxs = falses(size(midline_idxs))
    valid_idxs[minimum(find(midline_idxs))-1] = true
    valid_idxs[maximum(find(midline_idxs))+1] = true

    new_val = mean(v.data[valid_idxs, :, :], 1)
    for i in find(midline_idxs)
        v.data[i, :, :] = new_val
    end

    return v
end
