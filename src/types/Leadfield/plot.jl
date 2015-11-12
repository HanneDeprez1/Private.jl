using Winston

function plotSRF(X::AbstractVector, Y::AbstractVector, Z::AbstractVector; desiredPoints::Number=2048)

    # Reduce the number of points plotted so image size is reasonable
    newPoints = round(linspace(1,size(X)[1], desiredPoints))
    X = X[newPoints]
    Y = Y[newPoints]
    Z = Z[newPoints]

    p1 = FramedPlot(
        title="Back",
        xlabel="Left - Right",
        ylabel="Inferior - Superior")
    b1 = Points(X, Z)
    style(b1, kind="dot")
    add(p1, b1)

    p2 = FramedPlot(
        title="Side",
        xlabel="Posterior - Anterior",
        ylabel="Inferior - Superior")
    b2 = Points(Y, Z)
    style(b2, kind="dot")
    add(p2, b2)

    p3 = FramedPlot(
        title="Top",
        xlabel="Left - Right",
        ylabel="Posterior - Anterior")
    b3 = Points(X, Y)
    style(b3, kind="dot")
    add(p3, b3)

    t = Table(2,2)
    t[1,1] = p1
    t[1,2] = p2
    t[2,1] = p3

    return t

end
