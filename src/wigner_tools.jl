function wigner_normalise(q, p, W)
    """
    Normalise the Wigner function using the 2D integral
    """
    # Determine the norm
    N = int_2d(q, p, W)
    # Divide by norm
    return W ./ N
end

function calc_wigner_wqp(q, p, W)
    """
    Calculates the Wigner function in the q and p plane
    Wq(p) = int dq W(q, p)
    Wp(q) = int dp W(q, p)
    """

    # Find Wq
    Wq = [int_1d(p, W[i, :]) for i = 1:N]

    # Find Wp
    Wp = [int_1d(q, W[:, i]) for i = 1:N]

    return Wq, Wp
end

function calc_wigner_wqp_norm(q, p, W)
    """
    Calculates the norm of the Wigner function in the q and p plane
    """
    # Find Wq and Wp
    Wq, Wp = calc_wigner_wqp(q, p, W)

    # Calculate the norms
    q_norm = int_1d(q, Wq)
    p_norm = int_1d(p, Wp)

    return q_norm, p_norm
end

function calc_wigner_wqp_expect(q, p, W)
    """
    Calculates the expectation value of q and p
    <q> = int dq Wq(q) q
    <p> = int dp Wp(p) p
    """
    # Find W_q and W_p
    Wq, Wp = calc_wigner_wqp(q, p, W)

    # Calculate the expectation value
    q_ex = int_1d(q, Wq .* q)
    p_ex = int_1d(p, Wp .* p)

    return q_ex, p_ex
end


function calc_wigner_wqp2_expect(q, p, W)
    """
    Calculates the second moment of q and p
    <q^2> = int dq Wq(q) q^2
    <p^2> = int dp Wp(p) p^2
    """
    # Find Wq and Wp
    Wq, Wp = calc_wigner_wqp(q, p, W)

    # Calculate the second moment
    q2_ex = int_1d(q, Wq .* q .^ 2)
    p2_ex = int_1d(p, Wp .* p .^ 2)
    return q2_ex, p2_ex
end


function calc_wigner_uncertainty_principle(q, p, W)
    """
    Calculates the uncertainty principle
    sqrt(<q^2> - <q>^2) * sqrt(<p^2> - <p>^2) =? 1/2
    """
    # Find the q and p expectation values
    q_ex, p_ex = calc_wigner_wqp_expect(q, p, W)
    # Find the q and p second moments
    q2_ex, p2_ex = calc_wigner_wqp2_expect(q, p, W)
    return sqrt(q2_ex - q_ex^2) * sqrt(p2_ex - p_ex^2)
end

function calc_wigner_entropy(q, p, W)
    """
    !!!DOUBLE CHECK THIS!!!
    Calculates the entropy of the Wigner function
    S = -int dq dp W(q, p) log(W(q, p))
    """
    # Find the entropy
    return -int_2d(q, p, W .* log.(abs.(W)))
end