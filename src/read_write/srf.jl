"""
    readSRF(fname::AbstractString)

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
