
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
