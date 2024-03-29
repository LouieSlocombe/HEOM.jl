function wigner_normalise(q, p, W; k=5)
    """
    Normalise the Wigner function using the 2D integral
    """
    # Determine the norm
    norm = int_2d(q, p, W; k=k)
    # Divide by norm
    return W ./ norm
end

function calc_wigner_wqp(q, p, W; k=5)
    """
    Calculates the Wigner function in the q and p plane
    Wq(p) = ∫ W(q, p) dq
    Wp(q) = ∫ W(q, p) dp
    """
    # Get the number of q and p points
    Nq = length(q)
    Np = length(p)
    # Find Wq
    Wq = [int_1d(p, W[i, :]; k=k) for i = 1:Nq]
    # Find Wp
    Wp = [int_1d(q, W[:, i]; k=k) for i = 1:Np]
    return Wq, Wp
end

function calc_wigner_wqp_norm(q, p, W; k=5)
    """
    Calculates the norm of the Wigner function in the q and p plane
    """
    # Find Wq and Wp
    Wq, Wp = calc_wigner_wqp(q, p, W; k=k)
    # Calculate the norms
    q_norm = int_1d(q, Wq; k=k)
    p_norm = int_1d(p, Wp; k=k)
    return q_norm, p_norm
end

function calc_wigner_o_expect(q, p, W, O; k=5)
    """
    Calculates the expecation value of an operator O
    <O> = ∫∫ W(q, p) O(q, p) dq dp
    """
    # Integrate over q and p
    return int_2d(q, p, W .* O; k=k)
end

function calc_wigner_wqp_expect(q, p, W; k=5)
    """
    Calculates the expectation value of q and p
    <q> = ∫∫ W(q, p) q dq dp
    <p> = ∫∫ W(q, p) p dq dp
    """
    # Calculate the expectation value
    q_ex = calc_wigner_o_expect(q, p, W, q; k=k)
    p_ex = calc_wigner_o_expect(q, p, W, p; k=k)
    return q_ex, p_ex
end

function calc_wigner_wqp2_expect(q, p, W; k=5)
    """
    Calculates the second moment of q and p
    <q^2> = ∫∫ W(q, p) q^2 dq dp
    <p^2> = ∫∫ W(q, p) p^2 dq dp
    """
    # Calculate the second moment
    q2_ex = calc_wigner_o_expect(q, p, W, q .^ 2; k=k)
    p2_ex = calc_wigner_o_expect(q, p, W, p .^ 2; k=k)
    return q2_ex, p2_ex
end

function calc_wigner_uncertainty_principle(q, p, W; k=5)
    """
    Calculates the uncertainty principle
    sqrt(<q^2> - <q>^2) * sqrt(<p^2> - <p>^2) >= ħ/2
    """
    # Find the q and p expectation values
    q_ex, p_ex = calc_wigner_wqp_expect(q, p, W; k=k)
    # Find the q and p second moments
    q2_ex, p2_ex = calc_wigner_wqp2_expect(q, p, W; k=k)
    return sqrt(q2_ex - q_ex^2) * sqrt(p2_ex - p_ex^2)
end

function calc_wigner_entropy(q, p, W; k=5)
    """
    !!!DOUBLE CHECK THIS!!!
    Calculates the entropy of the Wigner function
    S = -∫∫ W(q, p) log(W(q, p)) dq dp
    """
    # Calculate the log of the Wigner function
    W_log = @. W * log(abs(W))
    # Find the entropy
    return -int_2d(q, p, W_log; k=k)
end

function calc_classical_hamiltonian(p, v, mass)
    """
    Calculates the classical Hamiltonian
    H = P^2 / (2m) + V
    """
    P = repeat(reshape(p, 1, :), length(p), 1)
    return @. P^2 / (2.0 * mass) + v
end

function calc_wigner_energy_expect(q, p, v, mass, W; k=5)
    """
    Calculates the energy expectation value
    <H> = ∫∫ W(q, p) H(q, p) dq dp
    """
    # Calculate the classical Hamiltonian
    H = calc_classical_hamiltonian(p, v, mass)
    # Calculate the expectation value
    return calc_wigner_o_expect(q, p, W, H; k=k)
end

function calc_wigner_purity(q, p, W; k=5)
    """
    Calculates the purity of the Wigner function
    P = 2πħ ∫∫ W(q, p)^2 dq dp
    """
    return 2.0 * pi * h_bar * calc_wigner_o_expect(q, p, W, W; k=k)
end

function calc_wigner_s2_entropy(q, p, W; k=5)
    """
    Calculates the entropy of the Wigner function
    S = 1 - 2πħ ∫∫ W(q, p)^2 dq dp
    """
    return 1.0 - calc_wigner_purity(q, p, W; k=k)
end

function calc_wigner_autocorrelation(q, p, sol; k=5)
    """
    Calculates the autocorrelation function
    C(t) = ∫∫ W(q, p, t) W(q, p, 0) dq dp
    """
    # Get the number of time steps
    nt = length(sol.t)
    # Get the Wigner function
    W = sol.u
    # Calculate the autocorrelation
    C = [int_2d(q, p, W[:, :, i] .* W[:, :, 1]; k=k) for i = 1:nt]
    # Normalise
    return C ./ C[1]
end

function calc_wigner_probability(q, p, q0, W; k=5, f_renorm=false)
    """
    Calculates the probability of finding the particle in a region of phase space
    P = ∫∫ W(q, p) h_q0(q, p) dq dp
    """
    if f_renorm
        W = wigner_normalise(q, p, W; k=k)
    end
    # Calculate the step function
    h_q0 = step_function(q, W, q0)
    # Calculate the probability
    return int_2d(q, p, h_q0; k=k)
end

function calc_wigner_k_qm(q, p, q0, sol; k=5, f_renorm=false)
    """
    Calculates the rate of change of the probability 
    of finding the particle in a region of phase space
    dP/dt = ∫∫ ∂W/∂t h_q0(q, p) dq dp
    """
    time = sol.t
    # Calculate the occupation probability
    dP_p = [calc_wigner_probability(q, p, q0, sol[i]; k=k, f_renorm=f_renorm) for i = 1:length(time)]
    # Calculate the rate of change
    dP_p_dt = derivative_1d_interp(time, dP_p, 1; k=k)
    # Prevent inf at zero
    #dP_p_dt[1] = dP_p_dt[2]
    return dP_p_dt
end