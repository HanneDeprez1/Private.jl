"""
    readELP(fname::AbstractString)

Import electrode data in the BESA elp format`

#### Arguments

* `fname` = location and name of file to import
* `verbose` = should the function output user feedback


name, a, b = readELP(verbose=false)
"""
function readELP(fname::AbstractString = Pkg.dir("Private", "test", "data", "Standard-10-10-Cap81.elp"); verbose::Bool=false)

    debug("Electrode file: $(fname)")

    df = readtable(fname, header = false, separator = ' ')

    if verbose
        println("Electrodes:  $(length(df[:x1]))")
    end

    return df[:x1].data, df[:x2].data, df[:x3].data
end
