function wigner_normalise(q, p, W)
    """
    Normalise the Wigner function using the 2D integral
    """
    # Determine the norm
    N = int_2d(q, p, W)
	# Divide by norm
    return W ./ N
end