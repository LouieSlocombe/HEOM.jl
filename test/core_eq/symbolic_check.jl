using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings
using Test
const H = HEOM

# Define the equation and its derivatives
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
# Check the the Wigner master equation
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=false)
# println("eq = ", eq)
# println(latexify(eq))
# println(H.latexify_nice(eq))
wm_eq = (-p * Differential(q)(W(q, p, t))) / m + Differential(q)(V(q, t)) * Differential(p)(W(q, p, t)) + (1 // 24) * (ħ^2) * Differential(q)(Differential(q)(Differential(q)(V(q, t)))) * Differential(p)(Differential(p)(Differential(p)(W(q, p, t))))
@test H.eq_simple(eq - wm_eq) == 0

