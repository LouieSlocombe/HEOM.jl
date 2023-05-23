# Symbolic form of the Wigner master equation
function generate_wm_eq(W, args; f_simple=false)
    # Kinetic term
    WM_kin = -args.p / args.m * args.Dq(W(args.q, args.p, args.t))
    # # First order potential term
    WM_pot = args.Dq(args.V(args.q, args.t)) * args.Dp(W(args.q, args.p, args.t))
    # # Third order potential term
    WM_pot_aux = -args.ħ^2 / 24 * args.Dqqq(args.V(args.q, args.t)) * args.Dppp(W(args.q, args.p, args.t))
    # # Make the quantum operator
    L_QM = WM_kin + WM_pot + WM_pot_aux
    if f_simple
        L_QM = eq_simple(L_QM)
    end
    return L_QM
end

# Symbolic form of the truncated Wigner master equation
function generate_wm_trunc_eq(W, args; f_simple=false)
    # Kinetic term
    WM_kin = -args.p / args.m * args.Dq(W(args.q, args.p, args.t))
    # First order potential term
    WM_pot = args.Dq(args.V(args.q, args.t)) * args.Dp(W(args.q, args.p, args.t))
    # Make the quantum operator
    L_QM = WM_kin + WM_pot
    if f_simple
        L_QM = eq_simple(L_QM)
    end
    return L_QM
end

# Symbolic form of the LL_HT_M master equation
function generate_LL_HT_M_eq(W, args; f_simple=false)
    # Kinetic term
    WM_kin = -args.p / args.m * args.Dq(W(args.q, args.p, args.t))
    # First order potential term
    WM_pot = args.Dq(args.V(args.q, args.t)) * args.Dp(W(args.q, args.p, args.t))
    # Third order potential term
    WM_pot_aux = -args.ħ^2 / 24 * args.Dqqq(args.V(args.q, args.t)) * args.Dppp(W(args.q, args.p, args.t))
    # Dissipation term
    Diss = args.ζ * args.Dp(args.p * W(args.q, args.p, args.t))
    # Decoherence term
    Deco = args.ζ * args.m / args.β * args.Dpp(W(args.q, args.p, args.t))
    # Make the quantum operator
    L_QM = WM_kin + WM_pot + WM_pot_aux + Diss + Deco
    if f_simple
        L_QM = eq_simple(L_QM)
    end
    return L_QM
end

# Symbolic form of the truncated LL_HT_M master equation
function generate_LL_HT_M_trunc_eq(W, args; f_simple=false)
    # Kinetic term
    WM_kin = -args.p / args.m * args.Dq(W(args.q, args.p, args.t))
    # First order potential term
    WM_pot = args.Dq(args.V(args.q, args.t)) * args.Dp(W(args.q, args.p, args.t))
    # Dissipation term
    Diss = args.ζ * args.Dp(args.p * W(args.q, args.p, args.t))
    # Decoherence term
    Deco = args.ζ * args.m / args.β * args.Dpp(W(args.q, args.p, args.t))
    # Make the quantum operator
    L_QM = WM_kin + WM_pot + Diss + Deco
    if f_simple
        L_QM = eq_simple(L_QM)
    end
    return L_QM
end

# Symbolic form of the LL_HT_NM master equation
function generate_LL_HT_NM_n_hierarchy(n, N, W, args; f_simple=false)
    # Makes the nth hierarchy
    i_n = n - 1
    #LW = generate_wm_eq(n)
    if n == 1
        # Chop off the terms below the zeroth hierarchy 
        hier = -(args.L + i_n * args.γ) * W(args.q, args.p, args.t)[n] + args.Dp(W(args.q, args.p, args.t)[n+1])
    elseif n == N
        # Chop off the terms above the Nth hierarchy
        hier = -(args.L + i_n * args.γ) * W(args.q, args.p, args.t)[n] + i_n * args.γ * args.ζ * args.p * W(args.q, args.p, args.t)[n-1] + i_n * args.γ * args.ζ * args.m / args.β * args.Dp(W(args.q, args.p, args.t)[n-1])
    else
        # All terms above and below
        hier = -(args.L + i_n * args.γ) * W(args.q, args.p, args.t)[n] + args.Dp(W(args.q, args.p, args.t)[n+1]) + i_n * args.γ * args.ζ * args.p * W(args.q, args.p, args.t)[n-1] + i_n * args.γ * args.ζ * args.m / args.β * args.Dp(W(args.q, args.p, args.t)[n-1])
    end
    # Check to simplify
    if f_simple
        hier = eq_simple(hier)
    end
    return hier
end