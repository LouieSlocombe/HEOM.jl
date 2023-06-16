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