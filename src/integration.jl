function int_1d(x, y; x0=nothing, x1=nothing, k=3)
    """
    Perform numerical integration by doing a 1D spline then integrate
    using the Dierckx package
    """
    # Find the domain to calculate the integral
    if x0 == nothing
        x0 = minimum(x)
    end
    if x1 == nothing
        x1 = maximum(x)
    end

    # Interpolate
    spl = Spline1D(x, y; k=k)
    # Intergrate
    return integrate(spl, x0, x1)
end

function int_2d(
    x,
    y,
    z;
    x0=nothing,
    x1=nothing,
    y0=nothing,
    y1=nothing,
    k=3,
)
    """
    Perform numerical integration by doing a 2D spline then integrate
    using the Dierckx package
    """
    # Find the domain to calculate the integral
    if x0 == nothing
        x0 = minimum(x)
    end
    if x1 == nothing
        x1 = maximum(x)
    end
    if y0 == nothing
        y0 = minimum(y)
    end
    if y1 == nothing
        y1 = maximum(y)
    end

    # Interpolate
    spl = Spline2D(x, y, z; kx=k, ky=k)
    # Intergrate
    return integrate(spl, x0, x1, y0, y1)
end