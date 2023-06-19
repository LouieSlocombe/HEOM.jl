function best_time_units(time)
    prefix_name = ["mili-", "micro-", "nano-", "pico-", "femto-", "atto-"]
    prefix_unit = ["m", L"\mu", "n", "p", "f", "a"]
    prefix_amount = [1.0e-3, 1.0e-6, 1.0e-9, 1.0e-12, 1.0e-15, 1.0e-18]

    # Convert to array and convert time from AUT to SI
    time = Array(time) .* au_time

    # Compare
    max_time = maximum(time)
    dif = @. abs(log10(max_time) - log10(prefix_amount))

    # Pick best range
    min_loc = argmin(dif)

    prefix = prefix_unit[min_loc]
    name = prefix_name[min_loc]
    # Convert the new time
    new_time = time ./ prefix_amount[min_loc]

    return new_time, prefix, name
end

function calculate_fwhm(x, y)
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
        left_margin=2Plots.mm
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