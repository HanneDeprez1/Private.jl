
using ProgressMeter


##########################
#
# Types - SSR
#
##########################


"""
Linearly constrained minimum variance (LCMV) beamformer for epoched data.
NAI is ratio between stimulus and control data.


### Input

* s = stimulus condition SSR data with epochs pre calculated
* n = control condition SSR with epochs pre calculated
* l = leadfield information

* foi = frequency of interest for cross power spectral density calculations
* fs = sample rate
* n_epochs = number of epochs to average down to with aim of reducing noise

"""
function beamformer_lcmv(s::SSR, n::SSR, l::Leadfield;
                         foi::Real=modulationrate(s), fs::Real=samplingrate(s), n_epochs::Int=0, kwargs...)

    Logging.info("Performing LCMV beamforming on signal with noise data as reference")

    if !haskey(s.processing, "epochs") || !haskey(n.processing, "epochs")
        Logging.critical("Epochs not calculated")
    end

    if n_epochs > 0
        s.processing["epochs"] = reduce_epochs(s.processing["epochs"], n_epochs)
        n.processing["epochs"] = reduce_epochs(n.processing["epochs"], n_epochs)
    end

    l = match_leadfield(l, s)

    V, N, NAI = beamformer_lcmv(s.processing["epochs"], n.processing["epochs"], l.L, fs, foi; kwargs...)

    VolumeImage(vec(NAI), "NAI", l.x, l.y, l.z, [1.0], "LCMV", Dict(), "Talairach")
end


##########################
#
# Low level functions
#
##########################



"""
Linearly constrained minimum variance (LCMV) beamformer for epoched data

LCMV beamformer returning neural activity index (NAI), source and noise variance (Van Veen et al 1997).
Source space projection is implemented (Sekihara et al 2001).


### Literature

Localization of brain electrical activity via linearly constrained minimum variance spatial filtering
Van Veen, B. D., van Drongelen, W., Yuchtman, M., & Suzuki, A.
Biomedical Engineering, IEEE Transactions on, 44(9):867–880, 1997.

Reconstructing spatio-temporal activities of neural sources using an meg vector beamformer technique
Kensuke Sekihara, Srikantan S Nagarajan, David Poeppel, Alec Marantz, and Yasushi Miyashita.
Biomedical Engineering, IEEE Transactions on, 48(7):760–771, 2001.


### Input

* x = M * T x N matrix = Signal M sample measurements over T trials on N electrodes
* n = M * T x N matrix = Noise  M sample measurements over T trials on N electrodes
* H = L x D x N matrix = Forward head model for L locations on N electrodes in D dimensions
* fs = Sample rate
* foi = Frequency of interest for cross spectral density
* freq_pm = Frequency above and below `foi` to include in csd calculation (1.0)

"""
function beamformer_lcmv{A <: AbstractFloat}(x::Array{A, 3}, n::Array{A, 3}, H::Array{A, 3}, fs::Real, foi::Real;
                         freq_pm::Real = 0.5, kwargs...)

    Logging.debug("Starting LCMV beamforming on epoch data of size $(size(x, 1)) x $(size(x, 2)) x $(size(x, 3)) and $(size(n, 1)) x $(size(n, 2)) x $(size(n, 3))")

    # Constants
    M = size(x, 1)   # Samples
    N = size(x, 3)   # Sensors
    L = size(H, 1)   # Locations
    D = size(H, 2)   # Dimensions

    # Check input
    @assert size(n, 3) == N     # Ensure inputs match
    @assert size(H, 3) == N     # Ensure inputs match
    @assert M > N               # Should have more samples than sensors
    @assert !any(isnan(x))
    @assert !any(isnan(n))
    @assert !any(isnan(H))

    Logging.debug("LCMV epoch beamformer using $M samples on $N sensors for $L sources over $D dimensions")

    C = cross_spectral_density(x, foi - freq_pm, foi + freq_pm, fs)
    Q = cross_spectral_density(n, foi - freq_pm, foi + freq_pm, fs)

    # More (probably unnecessary) checks
    @assert size(C) == size(Q)
    @assert C != Q

    beamformer_lcmv(C, Q, H; kwargs...)
end


function beamformer_lcmv{A <: AbstractFloat}(C::Array{Complex{A}, 2}, Q::Array{Complex{A}, 2}, H::Array{A, 3};
                              subspace::A=0.95, regularisation::A=0.003, progress::Bool=false, kwargs...)

    Logging.debug("Computing LCMV beamformer from CPSD data")

    N = size(C, 1)   # Sensors
    L = size(H, 1)   # Locations

    # Space to save results
    Variance  = Array(Float64, (L, 1))         # Variance
    Noise     = Array(Float64, (L, 1))         # Noise
    NAI       = Array(Float64, (L, 1))         # Neural Activity Index
    Logging.debug("Result variables pre allocated")

    # TODO before or after subspace?
    # Default as suggested in discussion of Sekihara
    if regularisation > 0
        S = svdfact(real(C)).S[1]
        C = C + regularisation * S * eye(C)
        Logging.debug("Regularised signal matrix with lambda = $(S * regularisation)")
        S = svdfact(real(Q)).S[1]
        Q = Q + regularisation * S * eye(Q)
        Logging.debug("Regularised noise matrix with lambda = $(S * regularisation)")
    end

    if subspace > 0

        # Create subspace from singular vectors
        ss, k = retain_svd(real(C), subspace)
        ss = ss'

        Logging.debug("Subspace constructed of $(size(ss, 1)) components constituting $(100*round(k, 5))% of power")

        # Apply subspace to signal and noise
        C = ss * C * ss'
        Q = ss * Q * ss'

        Logging.debug("Subspace projection calculated")

    else
        ss = eye(real(C))
    end

    # Compute inverse outside loop
    invC = pinv(C)
    invQ = pinv(Q)

    Logging.debug("Beamformer scan started")
    if progress; prog = Progress(L, 1, "  LCMV scan... ", 40); end
    for l = 1:L

        H_l = ss * squeeze(H[l,:,:], 1)'

        Variance[l], Noise[l], NAI[l] = beamformer_lcmv(invC, invQ, H_l)

        if progress; next!(prog); end
    end

    Logging.debug("Beamformer scan completed")

    return Variance, Noise, NAI
end


function beamformer_lcmv(invC::Array{Complex{Float64}, 2}, invQ::Array{Complex{Float64}, 2}, H::Array{Float64, 2})

    V_q = trace(pinv(H' * invC * H)[1:3, 1:3])   # Strength of source     Eqn 24: trace(3x3)

    N_q = trace(pinv(H' * invQ * H)[1:3, 1:3])   # Noise strength         Eqn 26: trace(3x3)

    NAI = V_q / N_q                              # Neural activity index  Eqn 27

    return abs(V_q), abs(N_q), abs(NAI)
end


###############################
#
# Cross power spectral density
#
###############################

"""
Compute complex cross spectral density of epoch data

Cross spectral density averaged over frequencies of interest

### Input

* epochs: data shaped as epochs (samples x trials x channels)
* fmin: minimum frequency of interest
* max: maximum frequency of interest

### Implementation

Currently uses MNE python library.
Will change to Synchrony.jl when its stabilised.

"""
function cross_spectral_density{T <: AbstractFloat}(epochs::Array{T, 3}, fmin::Real, fmax::Real, fs::Real; fsum::Bool=false, ignore::Real=Inf, ignore_pm::Real=1.0, kwargs...)

    @pyimport mne as mne
    @pyimport mne.time_frequency as tf

    Logging.debug("Cross power spectral density between $fmin - $fmax Hz for signal of $fs Hz")

    # Convert from EEG.jl to MNE format for epochs
    epochs = permutedims(epochs, [2, 3, 1])                      # Change to trials x channels x samples
    names = AbstractString[string(i) for i = 1:size(epochs, 2)]  # Hack to put in fake names
    events = ones(Int, size(epochs, 1), 3)                       # Make all events the same to use everything
    events[:, 1] = 0:size(events, 1) - 1

    # Run MNE processing
    i = pycall(mne.create_info, PyObject, ch_names = vec(names), ch_types = vec(repmat(["eeg"], size(epochs, 2))), sfreq = fs)
    epochs = mne.EpochsArray(epochs, i, events, 0.0, verbose = false)

    csd = tf.compute_epochs_csd(epochs, fmin = fmin, fmax = fmax, verbose = false, fsum = fsum)

    keep_idx = find(abs(AbstractFloat[c[:frequencies][1] for c in csd] - ignore) .> ignore_pm)
    csd = csd[keep_idx];
    Logging.debug("Averaging CSD over frequencies $(AbstractFloat[c[:frequencies][1] for c in csd])")
    a = zeros(csd[1][:data])
    for i in 1:length(csd)
        a = a .+ csd[i][:data]
    end
    a = a ./ length(csd)

    return a
end

function cross_spectral_density(s::SSR; fmin::Real=modulationrate(s)-0.5, fmax::Real=modulationrate(s)+0.5,
        fs::Real=samplingrate(s), kwargs...)
    cross_spectral_density(s.processing["epochs"], fmin, fmax, fs; kwargs...)
end


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





###########
#
# NEW
#
###########


"""
Return leadfield for target location
"""
function calculate_specific_leadfield2(ldf::Leadfield, target::Coordinate; kwargs...)

    calculate_specific_leadfield2(ldf, find_location(ldf, target); kwargs...)
end


"""
Return leadfield for target location index
"""
function calculate_specific_leadfield2(ldf::Leadfield, target::Int; kwargs...)

    calculate_specific_leadfield2(ldf.L, target; kwargs...)
end


"""
Return leadfield for target location with null placed at supress location
"""
function calculate_specific_leadfield2(ldf::Leadfield, target::Union{Coordinate, Int}, supress::Union{Coordinate, Int}; bilateral::Bool=false, kwargs...)

    if bilateral & isa(supress, Coordinate)
        supress.x = -1 * supress.x
    elseif bilateral & isa(supress, Int)
        supress = find_location(ldf, Talairach(-1 * ldf.x[supress], ldf.y[supress], ldf.z[supress]))
    end

    target_ldf = calculate_specific_leadfield2(ldf, target; kwargs...)
    supress_ldf = calculate_specific_leadfield2(ldf, supress; kwargs...)

    hcat(target_ldf, supress_ldf)
end


"""
Return leadfield for target location with nulls placed in supress region
"""
function calculate_specific_leadfield2(ldf::Leadfield, target::Union{Coordinate, Int}, supress::Union{Coordinate, Int}, region_size::Real; bilateral::Bool=false, kwargs...)

    if bilateral & isa(supress, Coordinate)
        supress.x = -1 * supress.x
    elseif bilateral & isa(supress, Int)
        supress = find_location(ldf, Talairach(-1 * ldf.x[supress], ldf.y[supress], ldf.z[supress]))
    end

    target_ldf = calculate_specific_leadfield2(ldf, target; kwargs...)
    supress_ldf = calculate_specific_leadfield_region(ldf, supress, region_size; kwargs...)

    hcat(target_ldf, supress_ldf )
end

function calculate_specific_leadfield2(ldf::Leadfield, target::Union{Coordinate, Int}, supress::Array{EEG.Talairach,1}, region_size::Real;
    bilateral::Bool=false, keep_vecs::Int=3, kwargs...)

    target_ldf = calculate_specific_leadfield2(ldf, target; kwargs...)

    for s in supress

        supress_ldf = calculate_specific_leadfield_region(ldf, s, region_size; kwargs...)

        supress_ldf = reduce_leadfield(supress_ldf, keep_vecs, keep_dims = 0; kwargs...)

        target_ldf = hcat(target_ldf, supress_ldf)
    end

    return target_ldf
end


"""
Return leadfield for location at given index
"""
function calculate_specific_leadfield2{A <: AbstractFloat}(ldf::Array{A, 3}, target::Int; kwargs...)

    squeeze(ldf[target,:,:], 1)'
end



"""
Return leadfield for given region given index
"""
function calculate_specific_leadfield_region(ldf::Leadfield, target::Int, region::Real; bilateral::Bool=false, kwargs...)

    supress_ldf = calculate_specific_leadfield_region(ldf, Talairach(ldf.x[target], ldf.y[target], ldf.z[target]), region; kwargs...)
end

"""
Return leadfield for given region given coordinaates
"""
function calculate_specific_leadfield_region(ldf::Leadfield, target::Coordinate, region::Real; bilateral::Bool=false, kwargs...)

    to_supress = falses(size(ldf.x))

    if bilateral
        target.x = -1 * target.x
    end

    for loc in 1:length(ldf.x)

        if euclidean(target, [ldf.x[loc], ldf.y[loc], ldf.z[loc]]) < region
            to_supress[loc] = true
        end
    end
    Ls = ldf.L[to_supress, :, :]
    Ls = reshape(permutedims(Ls, [3, 2, 1]), size(ldf.L, 3), size(Ls, 1) * 3)
end


"""
Reduce leadfield dimensionality
"""
function reduce_leadfield{F <: AbstractFloat}(ldf::Array{F, 2}, keep_vecs; keep_dims = 3, kwargs...)
    hcat(ldf[:, 1:keep_dims], svdfact(ldf[:, keep_dims+1:end]).U[:, 1:keep_vecs])

end






###########################
#
# Removed
#
###########################


"""
Linearly constrained minimum variance (LCMV) beamformer for epoched data.
NAI is ratio between frequeny of interest (`foi`) and specified control (`foi + noise_delta`) frequency.


### Input

* s = stimulus condition SSR data with epochs pre calculated
* l = leadfield information

* foi = frequency of interest for cross power spectral density calculations
* fs = sample rate
* n_epochs = number of epochs to average down to with aim of reducing noise

"""
function beamformer_lcmv(s::SSR, l::Leadfield;
                         foi::Real=modulationrate(s), fs::Real=samplingrate(s), n_epochs::Int=0,
                         freq_pm::Real = 0.5, noise_delta::Real=4.0, bilateral::Real = 15, kwargs...)

    Logging.info("Performing LCMV beamforming on signal data with $(foi + noise_delta) Hz as reference")

    if !haskey(s.processing, "epochs")
          Logging.critical("Epochs not calculated")
    end

    if n_epochs > 0
        s.processing["epochs"] = reduce_epochs(s.processing["epochs"], n_epochs)
    end

    l = match_leadfield(l, s)

    C = cross_spectral_density(s.processing["epochs"], foi - freq_pm, foi + freq_pm, fs; kwargs...)
    Q = cross_spectral_density(s.processing["epochs"], foi + noise_delta - freq_pm, foi + noise_delta + freq_pm, fs, ignore = foi; kwargs...)

    @assert size(C) == size(Q)
    @assert C != Q

    Logging.debug("Covariance matrices calculated and of size $(size(Q))")

    V, N, NAI = beamformer_lcmv(C, Q, l.L, l.x, l.y, l.z, bilateral; kwargs...)

    VolumeImage(vec(NAI), "NAI", l.x, l.y, l.z, [1.0], "LCMV", Dict(), "Talairach")
end
