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
    Calculates the expectation value of q and p
    <q> = int dq Wq(q) q
    <p> = int dp Wp(p) p
    """

    # Find W_q
    W_q = [int_1d(p, W[i, :]) for i = 1:N]

    # Find W_p
    W_p = [int_1d(q, W[:, i]) for i = 1:N]

    # Calculate the norms
    Q_norm = int_1d(q, W_q)
    P_norm = int_1d(p, W_p)

    # Calculate the expectation value
    Q_expect = int_1d(q, W_q .* q)
    P_expect = int_1d(p, W_p .* p)

    return W_q, W_p, Q_expect, P_expect, Q_norm, P_norm
end

function calc_wigner_wqp_expect(q, p, W)
    """
    Calculates the expecation value and the second moment for both the
    position and momentum
    <q> = int dq Wq(q) q
    <p> = int dp Wp(p) p
    <q^2> = int dq Wq(q) q^2
    <p^2> = int dp Wp(p) p^2
    """
    W_q, W_p, Q_expect, P_expect, _, _ = calc_wigner_wqp(q, p, W)
    Q2_expect = int_1d(q, W_q .* q .^ 2)
    P2_expect = int_1d(p, W_p .* p .^ 2)
    return Q2_expect, Q_expect, P2_expect, P_expect
end