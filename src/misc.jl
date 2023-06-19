function calculate_fwhm(x,y)
    """
    Calculate the full width at half maximum of a curve
    """
    y_hmax = maximum(y) * 0.5
    return x[findlast(y .>= y_hmax)] - x[findfirst(y .>= y_hmax)]
end

function step_function(x, y, x0)
    """
    Returns a step function with a step at x0
    """
    return @. 0.5 * (sign(x - x0) + 1.0) * y
end

function plot_general(x, y, xlab, ylab)
    """
    Plot general x y
    """
    p = plot(
        x,
        y,
        lw=2,
        legend=false,
        linecolor=:black,
        xlabel=xlab,
        ylabel=ylab,
    )
    return p
end

function linspace(x1, x2, n)
    """
    Julia version of np.linspace
    """
    return Array(collect(range(x1, x2, length=Int(n))))
end

function logspace(x1, x2, n)
    """
    Julia version of np.logspace
    Solution from
    https://stackoverflow.com/questions/59985035/does-there-exist-any-alternative-of-logspace-in-julia-v1-3-1
    """
    return Array(
        collect((10^y for y in range(log10(x1), log10(x2), length=Int(n)))),
    )
end