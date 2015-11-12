function import_leadfield(;verbose::Bool=false)

    L                     = readLFT(verbose=verbose)
    positions, neighbours = readLOC(verbose=verbose)
    names, a, b           = readELP(verbose=verbose)

    L = convert(Array{Float64, 3}, L)

    if size(L, 1) != size(positions,2); error("Mismatched leadfield and locations"); end
    if size(L, 2) != size(positions,1); error("Mismatched leadfield and locations"); end
    if size(L, 3) != size(names,1); error("Mismatched leadfield and names"); end

    Leadfield(L, vec(positions[1,:]), vec(positions[2,:]), vec(positions[3,:]), vec(names))
end
