"""
    plotSRF(X::AbstractVector, Y::AbstractVector, Z::AbstractVector)

Very simple plotting of surface results
"""
function plotSRF(X::AbstractVector, Y::AbstractVector, Z::AbstractVector)

    # TODO: remove Plots. once all other plotting code is gone

    a = Plots.scatter(X, -Y)
    b = Plots.scatter(Z, -Y)
    c = Plots.scatter(Z, X)

    Plots.plot(a, b, c, layout = @layout([a b c]))

end
