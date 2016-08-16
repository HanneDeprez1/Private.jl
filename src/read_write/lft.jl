"""
    readLFT(fname::AbstractString)

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
