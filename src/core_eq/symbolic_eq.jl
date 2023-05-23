# Symbolic form of the Wigner master equation
function generate_wm_eq(W, args; f_simple=false)
    # Kinetic term
    WM_kin = -args.p / args.m * args.Dq(W(args.q, args.p, args.t))
    # # First order potential term
    WM_pot = args.Dq(args.V(args.q, args.t)) * args.Dp(W(args.q, args.p, args.t))
    # # Third order potential term
    WM_pot_aux = args.ħ^2 / 24 * args.Dqqq(args.V(args.q, args.t)) * args.Dppp(W(args.q, args.p, args.t))
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
    WM_kin = -args.p / args.m * args.Dq(W)
    # First order potential term
    WM_pot = args.Dq(args.V) * args.Dp(W)
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
    WM_kin = -args.p / args.m * args.Dq(W)
    # First order potential term
    WM_pot = args.Dq(args.V) * args.Dp(W)
    # Third order potential term
    WM_pot_aux = args.ħ^2 / 24 * args.Dqqq(args.V) * args.Dppp(W)
    # Dissipation term
    Diss = args.ζ * args.Dp(args.p * W)
    # Decoherence term
    Deco = args.ζ * args.m / args.β * args.Dpp(W)
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
    WM_kin = -args.p / args.m * args.Dq(W)
    # First order potential term
    WM_pot = args.Dq(args.V) * args.Dp(W)
    # Dissipation term
    Diss = args.ζ * args.Dp(args.p * W)
    # Decoherence term
    Deco = args.ζ * args.m / args.β * args.Dpp(W)
    # Make the quantum operator
    L_QM = WM_kin + WM_pot + Diss + Deco
    if f_simple
        L_QM = eq_simple(L_QM)
    end
    return L_QM
end