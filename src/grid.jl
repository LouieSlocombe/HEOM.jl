function create_basis_even(n, q_range, p_range)
    # Calculate the step size
    dq = q_range / n
    dp = p_range / n

    # Make vectors
    q_vec = range(start=-q_range / 2, stop=q_range / 2 - dq, length=n) |> Array
    p_vec = range(start=-p_range / 2, stop=p_range / 2 - dp, length=n) |> Array

    # Make grid
    Q = repeat(q_vec, 1, n)
    P = repeat(reshape(p_vec, 1, :), n, 1)

    return q_vec, p_vec, Q, P, dq, dp
end