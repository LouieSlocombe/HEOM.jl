using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings
using Test
const H = HEOM

# Define the paramters and derivatives
@parameters t q p γ m β ω ħ ζ σ q0 p0 L
@variables W(..) V(..)
Dt = Differential(t)
Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

####################################################################################
# Check the Wigner master equation
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=false)
# println("eq = ", eq)
# println(latexify(eq))
# println(H.latexify_nice(eq))
eq_check = (-p * Differential(q)(W(q, p, t))) / m + Differential(q)(V(q, t)) * Differential(p)(W(q, p, t)) - (1 // 24) * (ħ^2) * Differential(q)(Differential(q)(Differential(q)(V(q, t)))) * Differential(p)(Differential(p)(Differential(p)(W(q, p, t))))
@test H.eq_simple(eq - eq_check) == 0

#####################################################################################
# Check the formalism allows for a potential and cancels out terms
v(q,t) = 0.5 * m * ω^2 * q^2
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=v, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=true)
eq_check = (-p * Differential(q)(W(q, p, t))) / m + m * q * (ω^2) * Differential(p)(W(q, p, t))
println("eq = ", eq)
println(latexify(eq))
println(H.latexify_nice(eq))
@test H.eq_simple(eq - eq_check) == 0

####################################################################################
# Check the truncated Wigner master equation
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, V=V)
eq = H.generate_wm_trunc_eq(W, args; f_simple=false)
eq_check = (-p * Differential(q)(W(q, p, t))) / m + Differential(q)(V(q, t)) * Differential(p)(W(q, p, t))
@test H.eq_simple(eq - eq_check) == 0

####################################################################################
# Check the LL_HT_M master equation
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, ζ=ζ, β=β, V=V, Dpp=Dpp, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_LL_HT_M_eq(W, args; f_simple=false)
eq_check = ζ * Differential(p)(p * W(q, p, t)) + (-p * Differential(q)(W(q, p, t))) / m + (m * ζ * Differential(p)(Differential(p)(W(q, p, t)))) / β + Differential(q)(V(q, t)) * Differential(p)(W(q, p, t)) - (1 // 24) * (ħ^2) * Differential(q)(Differential(q)(Differential(q)(V(q, t)))) * Differential(p)(Differential(p)(Differential(p)(W(q, p, t))))
@test H.eq_simple(eq - eq_check) == 0

####################################################################################
# Check the truncated LL_HT_M master equation
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, ζ=ζ, β=β, V=V, Dpp=Dpp)
eq = H.generate_LL_HT_M_trunc_eq(W, args; f_simple=false)
eq_check = ζ * Differential(p)(p * W(q, p, t)) + (-p * Differential(q)(W(q, p, t))) / m + (m * ζ * Differential(p)(Differential(p)(W(q, p, t)))) / β + Differential(q)(V(q, t)) * Differential(p)(W(q, p, t))
@test H.eq_simple(eq - eq_check) == 0

####################################################################################
# Check the LL_HT_NM terms
N_rank = 3
@variables W(..)[1:N_rank]
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, ζ=ζ, γ=γ, β=β, V=V, Dpp=Dpp, L=L)
n = 1
eq = H.generate_LL_HT_NM_n_hierarchy(n, N_rank, W, args; f_simple=false)
eq_check = Differential(p)((W(q, p, t))[2]) - L * (W(q, p, t))[1]
@test H.eq_simple(eq - eq_check) == 0
n = 2
eq = H.generate_LL_HT_NM_n_hierarchy(n, N_rank, W, args; f_simple=false)
eq_check = (-L - γ) * (W(q, p, t))[2] + (m * γ * ζ * Differential(p)((W(q, p, t))[1])) / β + p * γ * ζ * (W(q, p, t))[1] + Differential(p)((W(q, p, t))[3])
@test H.eq_simple(eq - eq_check) == 0