using DataFrames

"""
Import a leadfield matrix in the BESA LFT format

#### Arguments

* `fname`   = location and name of file to import
* `verbose` = should the function output user feedback

#### Returns

* Array representing leadfield values

#### Example

L = readLFT(verbose=true)
"""
function readLFT(fname::AbstractString = Pkg.dir("Private", "test", "data", "FEM_Default50MRI-Aniso3_CR80.lft");
        verbose::Bool = false)

    debug("Leadfield file: $(fname)")

    # Open file
    fid = open(fname, "r")

    # Read header information from file
    version         = read(fid, Int32, 1)[1]
    numSensors      = read(fid, Int32, 1)[1] # n_s
    numNodes        = read(fid, Int32, 1)[1] # n_n
    numOrientations = read(fid, Int32, 1)[1] # n_o

    # Read in maximum leadfield values
    numMaxLfValues = numNodes * numOrientations
    Lmax_jk = read(fid, Float32, numMaxLfValues)

    # Read in leadfield fractions
    numLfFracs = numSensors * numNodes[1] * numOrientations
    Lc_ijk = read(fid, Int16, numLfFracs)

    # Check that everything has been read
    if ~eof(fid)
        println("!! Something went wrong. File is longer than expected")
    end

    close(fid)

    ##  Convert to true leadfield values from fraction

    # Reshape Lmax
    # reshape requires particular inputs. Hopefully this is relaxed to save conversion
    Lmax_jk = reshape(transpose(Lmax_jk), (
        convert(Int64,numNodes),
        convert(Int64,numOrientations)))

    # Reshape Lc
    Lc_ijk = reshape(transpose(Lc_ijk), (
        convert(Int64, numNodes),
        convert(Int64, numOrientations),
        convert(Int64, numSensors)))

    # Scale leadfield
    # Lc_ijk = L_ijk * (SHRT_MAX / Lmax_jk) -> L_ijk = Lc_ijk / (SHRT_MAX / Lmax_jk)
    #                                          L_ijk = Lc_ijk * (Lmax_jk / SHRT_MAX)
    SHRT_MAX = 2^15-1

    L_ijk = Array(Float32, (
        convert(Int64, numNodes),
        convert(Int64, numOrientations),
        convert(Int64, numSensors)
        ))

    i = 1
    while i <= numNodes

        j = 1
        while j <= numOrientations

            k = 1
            while k <= numSensors

                L_ijk[i,j,k] =  Lc_ijk[i,j,k] * (Lmax_jk[i,j] / SHRT_MAX)

                k = k + 1
            end
            j = j + 1
        end
        i = i + 1
    end

    if verbose
        println("LFT version:            $(version)")
        println("Number of sensors:      $(numSensors)")
        println("Number of nodes:        $(numNodes)")
        println("Number of orientations: $(numOrientations)")
        println("Number of maximums:     $numMaxLfValues")
        println("Number of fractions:    $numLfFracs")
        println("Lc_ijk  = [$(size(Lc_ijk)[1]), $(size(Lc_ijk)[2]), $(size(Lc_ijk)[3])]")
        println("Lmax_jk = [$(size(Lmax_jk)[1]), $(size(Lmax_jk)[2])]")
        println("L_ijk   = [$(size(L_ijk)[1]), $(size(L_ijk)[2]), $(size(L_ijk)[3])]")
    end

    return L_ijk
end


"""
Import a leadfield location matrix in the BESA loc format

#### Arguments

* `fname` = location and name of file to import
* `verbose` = should the function output user feedback


#### Returns

locations = readLOC(verbose=true)
"""
function readLOC(fname::AbstractString = Pkg.dir("Private", "test", "data", "FEM_Default50MRI.loc");
    verbose::Bool=false)

    debug("Location file: $(fname)")

    # Open file
    fid = open(fname, "r")

    # Read header information from file
    version         = read(fid, Int32, 1)[1]
    numNodes        = read(fid, Int32, 1)[1]

    # Read in positions
    numPositions = 3*numNodes
    positions = read(fid, Float64, numPositions)

    # Read in number of neighbours
    numNeighbours = read(fid, Int32, numNodes)

    maxNeighbours = maximum(numNeighbours)
    neighbours = zeros(Int32, (convert(Int64, numNodes), convert(Int64,maxNeighbours)))

    i = 1
    while i <= numNodes
        j = 1
        while j <= numNeighbours[i]
            neighbours[i,j]= read(fid, Int32, 1)[1]
            j = j + 1
        end
        i = i + 1
    end

    # Check that everything has been read
    if ~eof(fid)
        println("!! Something went wrong. File is longer than expected")
    end

    close(fid)

    # Reshape position data
    positions = reshape(vec(positions), (Int(3), Int(numNodes)))
    positions = positions*1e3           # Convert to mm
    positions = positions .- 127.5
    positions = reshape(positions, (3, Int(numNodes)))

    if verbose
        println("LFT version:            $(version)")
        println("Number of nodes:        $(numNodes)")
        println("Number of positions:    $(numPositions)")
        println("Maximum # neighbours:   $(maxNeighbours)")

    end

    return positions, neighbours
end


"""
Import surface data in the BESA loc format

fname   = location and name of file to import
verbose = should the function output user feedback

surface = readSRF(verbose=false)
"""
function readSRF(fname::AbstractString = Pkg.dir("Private", "test", "data", "Default50Brain.srf"); verbose::Bool=false)

    debug("Surface file: $(fname)")

    # Open file
    fid = open(fname, "r")

    # Read header information from file
    version         = read(fid, Float32, 1)[1]
    reserved        = read(fid, Int32, 1)[1]
    numVertices     = read(fid, Int32, 1)[1]
    numTriangles    = read(fid, Int32, 1)[1]
    MeshCenterX     = read(fid, Float32, 1)[1]
    MeshCenterY     = read(fid, Float32, 1)[1]
    MeshCenterZ     = read(fid, Float32, 1)[1]
    Xcords          = read(fid, Float32, numVertices)
    Ycords          = read(fid, Float32, numVertices)
    Zcords          = read(fid, Float32, numVertices)
    Xcomps          = read(fid, Float32, numVertices)
    Ycomps          = read(fid, Float32, numVertices)
    Zcomps          = read(fid, Float32, numVertices)
    rConvColor      = read(fid, Float32, 1)[1]
    gConvColor      = read(fid, Float32, 1)[1]
    bConvColor      = read(fid, Float32, 1)[1]
    aConvColor      = read(fid, Float32, 1)[1]
    rConcColor      = read(fid, Float32, 1)[1]
    gConcColor      = read(fid, Float32, 1)[1]
    bConcColor      = read(fid, Float32, 1)[1]
    aConcColor      = read(fid, Float32, 1)[1]
    meshColor       = read(fid, Int32, numVertices)

    if verbose
        println("Version:            $(version)")
        println("Reserved:           $(reserved)")
        println("Vertices #:         $(numVertices)")
        println("Triangles #:        $(numTriangles)")

    end

    i = 1
    while i <= numVertices
        number = read(fid, Int32, 1)[1]
        # Just overwriting. Not saving data. TODO: Fix this
        neighbour = read(fid, Int32, number)
        i = i + 1
    end

    Triangles = read(fid, Int32, 3*numTriangles)
    numTriStripElements = read(fid, Int32, 1)[1]
    seqStripElements = read(fid, Int32, numTriStripElements)

    close(fid)

    return Xcords, Ycords, Zcords

end

"""
Import electrode data in the BESA elp format`

#### Arguments

* `fname` = location and name of file to import
* `verbose` = should the function output user feedback


name, a, b = readELP(verbose=false)
"""
function readELP(fname::AbstractString = Pkg.dir("Private", "test", "data", "Standard-10-10-Cap81.elp"); kwargs...)

    debug("Electrode file: $(fname)")

    file = open(fname, "r")

    name = AbstractString[]
    a = AbstractFloat[]
    b = AbstractFloat[]
    regexp = r"(EEG|POL)\s+(\w+'?)(\s+(-?\d*.\d+)?\s+(-?\d*.\d+))?"

    while !eof(file)

        dm = match(regexp, readline(file))

        if dm.captures[5] != nothing

            push!(name, dm.captures[1])
            push!(name, dm.captures[2])
            push!(name, dm.captures[3])
        end

    end

    return name, a, b
end
