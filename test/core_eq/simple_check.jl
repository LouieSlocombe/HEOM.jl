using HEOM, Symbolics, ModelingToolkit, LinearAlgebra
using Test
const H = HEOM

# Set things up

# Make the grid
n = 2^9
q_range = 12.0
p_range = 60.0
q_vec, p_vec, Q, P, dq, dp = H.create_basis_even(n, q_range, p_range)

# Define the equation and its derivatives
@parameters t q p q0 p0
@variables W(..)

# Define 2d equation
eq2d = exp(-q^2 - 0.05 * p^2)
eq2d_vec = H.make_discretised_2d(eq2d, [q, p], q_vec, p_vec, [])

####################################################################################
# Check the equation works for a simple case