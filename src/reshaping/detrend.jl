
"""
    savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

Polynomial smoothing with the Savitsky Golay filters.

Copied from https://github.com/blakejohnson/Qlab.jl/blob/master/src/SavitskyGolay.jl

#### Sources

Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf  
Python example: http://wiki.scipy.org/Cookbook/SavitzkyGolay
"""
function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

    #Some error checking
    @assert isodd(windowSize) "Window size must be an odd integer."
    @assert polyOrder < windowSize "Polynomial order must me less than window size."

    halfWindow = round(Int, (windowSize-1) / 2)

    #Setup the S matrix of basis vectors.
    S = zeros(windowSize, polyOrder+1)
    for ct = 0:polyOrder
        S[:,ct+1] = collect(-halfWindow:halfWindow).^(ct)
    end

    #Compute the filter coefficients for all orders
    #From the scipy code it seems pinv(S) and taking rows should be enough
    G = S*pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv);

    #Pad the signal with the endpoints and convolve with filter
    paddedX = [x[1]*ones(halfWindow); x; x[end]*ones(halfWindow)]
    y = conv(filterCoeffs[end:-1:1], paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]

end

"""
    trend(::SSR, channel::Int)

Determine the trend for a steady state response signal.
"""
function trend(a::SSR, chan::Int; secs = 0.5, order = 2)

    ls = size(a.data, 1)

    windowSize = round(Int, 8192 * secs)
    if iseven(windowSize)
        windowSize += 1
    end

    savitsky_golay(vec(a.data[:, chan]), windowSize, order)
end

"""
    detrend(::SSR)

Remove the trend from a steady state response.
"""
function detrend(a::SSR; kwargs...)

    for c in 1:size(a.data, 2)
        a.data[:, c] = a.data[:, c] - trend(a, c; kwargs...)
    end

    return a
end
