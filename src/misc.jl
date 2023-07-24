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

function printout(x)
    """
    Flushes out print statements so that they are seen on clusters
    """
    println(x)
    flush(stdout)
end

function calc_classical_rate(dG, T)
    """
    Calculates the classical rate
    k_cl = k_b T / h_bar * exp(-dG / (k_b T))
    dG: Barrier height in Hartree
    temperature: Temperature in Kelvin
    """
    # Si k_b used here so that the [s] units is in SI not AUT
    pre_term = @. (si_k_b * T) / (2.0 * pi * si_h_bar)
    # Put it all together, Hartree k_b used here to cancel with E (Hartree)
    return @. pre_term * exp(-dG / (k_b * T))
end

function findlocalmaxima(signal::Vector)
    """
    Taken from:
    https://discourse.julialang.org/t/how-to-identify-local-maxima-peaks-in-a-time-signal/6000
    Note this can be duped if there are two values close to each other which are minima
    """
    inds = Int[]
    if length(signal) > 1
        if signal[1] > signal[2]
            push!(inds, 1)
        end
        for i = 2:length(signal)-1
            if signal[i-1] < signal[i] > signal[i+1]
                push!(inds, i)
            end
        end
        if signal[end] > signal[end-1]
            push!(inds, length(signal))
        end
    end
    return inds
end

function find_double_well_barrier_loc(pot)
    # Find the location of two minima
    vals = findlocalmaxima(-pot)
    if length(vals) == 1
        # Assume failed to find two minima so might just have a barrier
        max_val, max_loc = findmax(pot)
    else
        # Find the maximum location of the peak between the two minima
        max_val, max_loc = findmax(pot[minimum(vals):maximum(vals)])
        # Correct the index
        max_loc = max_loc .+ minimum(vals)
    end
    return max_loc
end


function make_small_h(q_vec, v)
    # Make heaviside step function
    T = eltype(q_vec)
    # Find the barrier grid location
    max_loc = find_double_well_barrier_loc(v)

    Nq = length(q_vec)
    h = zeros(T, Nq)
    # Loop over the indexes
    for i = 1:Nq
        if i >= max_loc # if it is in the product region
            h[i] = 1.0 # set to one
        end
    end
    return h
end

function calc_local_spring(mass, v, q_vec)
    # Find the second derivative
    d2v_dq2 = derivative_1d_interp(q_vec, v, 2; k=5)
    # Get omega
    omega = @. sqrt(abs(d2v_dq2 / mass))
    return omega
end

function determine_best_local_spring(mass, v, q_vec)
    # Find the location of global minimum
    _, min_loc = findmin(v)
    # Get the spring constant
    omega = calc_local_spring(mass, v, q_vec)
    return omega[min_loc]
end

function rnd(x; n=3)
    """
    Rounds a number to n digits (3)
    """
    return round(x; sigdigits=n)
end

function trunc_pot(x; a=0.1)
    return @. a * atan(x / a)
end