using HEOM, Symbolics, ModelingToolkit, LinearAlgebra
using Test
const H = HEOM

tol = 1e-2

# Make the grid
n = 2^9
q_range = 12.0
p_range = 60.0
q_vec, p_vec, Q, P, dq, dp = H.create_basis_even(n, q_range, p_range)

# Define the equation and its derivatives
@parameters t q p q0 p0
@variables W(..)

eq = exp(-q^2)
eq_vec = H.make_discretised_1d(eq, [q], q_vec, [])

Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

# Define 2d equation
eq2d = exp(-q^2 - 0.05 * p^2)
eq2d_vec = H.make_discretised_2d(eq2d, [q, p], q_vec, p_vec, [])

# 1D first order derivative
d1_eq1d = expand_derivatives(Dq(eq))
d1_eq1d_vec = H.make_discretised_1d(d1_eq1d, [q], q_vec, [])

####################################################################################
# 1D first order derivative using FFT
Dq_fft = H.prepare_fft_der(q_vec, q_range; order=1)
# Older version
fft_der = zeros(n)
@test ≈(Dq_fft * eq_vec, d1_eq1d_vec; atol=tol)
mul!(fft_der, Dq_fft, eq_vec)
@test ≈(fft_der, d1_eq1d_vec; atol=tol)

# 1D first order derivative using finite difference
d_dq = H.prepare_fd_dq(1, 2, dq, n)
fd_der = zeros(n)
H.dq_fd!(fd_der, eq_vec, d_dq)
@test ≈(fd_der, d1_eq1d_vec; atol=tol)

d_dq = H.prepare_fd_dq(1, 4, dq, n)
fd_der = zeros(n)
H.dq_fd!(fd_der, eq_vec, d_dq)
@test ≈(fd_der, d1_eq1d_vec; atol=tol)

####################################################################################
# 1D second order derivative
d2_eq1d = expand_derivatives(Dqq(eq))
d2_eq1d_vec = H.make_discretised_1d(d2_eq1d, [q], q_vec, [])

Dqq_fft = H.prepare_fft_der(q_vec, q_range; order=2)
@test ≈(Dqq_fft * eq_vec, d2_eq1d_vec; atol=tol)
fft_der = zeros(n)
mul!(fft_der, Dqq_fft, eq_vec)
@test ≈(fft_der, d2_eq1d_vec; atol=tol)

d_dqq = H.prepare_fd_dq(2, 2, dq, n)
fd_der = zeros(n)
H.dq_fd!(fd_der, eq_vec, d_dqq)
@test ≈(fd_der, d2_eq1d_vec; atol=tol)

####################################################################################
# 2D first order derivative d/dq
d1_eq2d = expand_derivatives(Dq(eq2d))
d1_eq2d_vec = H.make_discretised_2d(d1_eq2d, [q, p], q_vec, p_vec, [])

ft, FDq, ik = H.prepare_fft_dq(q_vec, eq2d_vec; order=1)
fft_der = zeros(n, n)
H.dq_fft!(fft_der, eq2d_vec, ft, FDq, ik)
@test ≈(fft_der, d1_eq2d_vec; atol=tol)

d_dq = H.prepare_fd_dq(1, 4, dq, n)
fd_der = zeros(n, n)
H.dq_fd!(fd_der, eq2d_vec, d_dq)
@test ≈(fd_der, d1_eq2d_vec; atol=tol)

####################################################################################
# 2D first order derivative d/dp
d1_eq2d = expand_derivatives(Dp(eq2d))
d1_eq2d_vec = H.make_discretised_2d(d1_eq2d, [q, p], q_vec, p_vec, [])

ft, FDp, ik = H.prepare_fft_dp(p_vec, eq2d_vec; order=1)
fft_der = zeros(n, n)
H.dp_fft!(fft_der, eq2d_vec, ft, FDp, ik)
@test ≈(fft_der, d1_eq2d_vec; atol=tol)

d_dp = H.prepare_fd_dp(1, 4, dp, n)
fd_der = zeros(n, n)
H.dp_fd!(fd_der, eq2d_vec, d_dp)
@test ≈(fd_der, d1_eq2d_vec; atol=tol)




####################################################################################
# 2D d2/dq2 + d2/dp2 + d/dq + d/dp
d2_eq2d = expand_derivatives(Dqq(eq2d) + Dpp(eq2d)+ Dq(eq2d) + Dp(eq2d))
d2_eq2d_vec = H.make_discretised_2d(d2_eq2d, [q, p], q_vec, p_vec, [])

ft, FDq, ik = H.prepare_fft_dq(q_vec, eq2d_vec; order=2)
d2_dq2 = zeros(n, n)
H.dq_fft!(d2_dq2, eq2d_vec, ft, FDq, ik)

ft, FDp, ik = H.prepare_fft_dp(p_vec, eq2d_vec; order=2)
d2_dp2 = zeros(n, n)
H.dp_fft!(d2_dp2, eq2d_vec, ft, FDp, ik)

ft, FDq, ik = H.prepare_fft_dq(q_vec, eq2d_vec; order=1)
d1_dq = zeros(n, n)
H.dq_fft!(d1_dq, eq2d_vec, ft, FDq, ik)
ft, FDq, ik = H.prepare_fft_dp(p_vec, eq2d_vec; order=1)
d1_dp = zeros(n, n)
H.dp_fft!(d1_dp, eq2d_vec, ft, FDp, ik)

fft_der = @. d2_dq2 + d2_dp2 + d1_dq + d1_dp
@test ≈(fft_der, d2_eq2d_vec; atol=tol)

d_dq = H.prepare_fd_dq(1, 4, dq, n)
d_dqq = H.prepare_fd_dq(2, 4, dq, n)
d_dp = H.prepare_fd_dp(1, 4, dp, n)
d_dpp = H.prepare_fd_dp(2, 4, dp, n)

dw_dq = zeros(n, n)
dw_dqq = zeros(n, n)
dw_dp = zeros(n, n)
dw_dpp = zeros(n, n)

H.dq_fd!(dw_dq, eq2d_vec, d_dq)
H.dq_fd!(dw_dqq, eq2d_vec, d_dqq)
H.dp_fd!(dw_dp, eq2d_vec, d_dp)
H.dp_fd!(dw_dpp, eq2d_vec, d_dpp)

fd_der = @. dw_dq + dw_dp + dw_dqq + dw_dpp
@test ≈(fd_der, d2_eq2d_vec; atol=tol)