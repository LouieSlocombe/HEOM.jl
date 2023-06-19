function w0_thermal(q, p, v, mass, beta)
    """
    Calculates the Wigner function for a thermal state
    """
    # Calculate the classical Hamiltonian
    ham = calc_classical_hamiltonian(p, v, mass)
    # Calculate the Wigner function
    W0 = @. exp(-beta * ham)
    # Normalise 
    return wigner_normalise(q, p, W0)
end

function w0_thermal_step(q, p, v, mass, beta, q0)
    """
    Calculates the Wigner function for a thermal state but with a step function
    """
    # Calculate the Wigner function
    W0 = w0_thermal(q, p, v, mass, beta)
    # Calculate the step function
    W0 = step_function(q, W0, q0)
    # Normalise
    return wigner_normalise(q, p, W0)
end

function w0_gaussian(q, p, q0, p0, σ, ħ)
    """
    Calculates the Wigner function for a Gaussian state
    """
    # Make grid
    n = length(q)
    Q = repeat(q, 1, n)
    P = repeat(reshape(p, 1, :), n, 1)
    # Calculate the Wigner function
    W0 = @. 1 / (π * ħ) * exp(-2 * σ * (Q - q0)^2 - 1 / (2 * ħ^2 * σ) * (P - p0)^2)
    # Normalise
    return wigner_normalise(q, p, W0)
end