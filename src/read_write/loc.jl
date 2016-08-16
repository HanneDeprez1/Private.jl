"""
    readLOC(fname::AbstractString)

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
