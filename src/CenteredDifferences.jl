module CenteredDifferences
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
using SciMLBase: AbstractDiffEqLinearOperator
import Base: +, -, *, /, \, size, getindex, setindex!, Matrix, convert, ==

# All of this is taken from https://github.com/SciML/DiffEqOperators.jl
# Which is under the MIT license
# All credit goes to the authors of that package!
# Go check it out!
# Copyright (c) 2017: shivin9.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

abstract type AbstractDiffEqAffineOperator{T} end
abstract type AbstractDerivativeOperator{T} <: AbstractDiffEqLinearOperator{T} end
abstract type AbstractMatrixFreeOperator{T} <: AbstractDiffEqLinearOperator{T} end

function calculate_weights(order::Int, x0::T, x::AbstractVector) where {T<:Real}
    #=
        order: The derivative order for which we need the coefficients
        x0   : The point in the array 'x' for which we need the coefficients
        x    : A dummy array with relative coordinates, e.g., central differences
               need coordinates centred at 0 while those at boundaries need
               coordinates starting from 0 to the end point
        The approximation order of the stencil is automatically determined from
        the number of requested stencil points.
    =#
    N = length(x)
    @assert order < N "Not enough points for the requested order."
    M = order
    c1 = one(T)
    c4 = x[1] - x0
    C = zeros(T, N, M + 1)
    C[1, 1] = 1
    @inbounds for i in 1:(N-1)
        i1 = i + 1
        mn = min(i, M)
        c2 = one(T)
        c5 = c4
        c4 = x[i1] - x0
        for j in 0:(i-1)
            j1 = j + 1
            c3 = x[i1] - x[j1]
            c2 *= c3
            if j == i - 1
                for s in mn:-1:1
                    s1 = s + 1
                    C[i1, s1] = c1 * (s * C[i, s] - c5 * C[i, s1]) / c2
                end
                C[i1, 1] = -c1 * c5 * C[i, 1] / c2
            end
            for s in mn:-1:1
                s1 = s + 1
                C[j1, s1] = (c4 * C[j1, s1] - s * C[j1, s]) / c3
            end
            C[j1, 1] = c4 * C[j1, 1] / c3
        end
        c1 = c2
    end
    #=
        This is to fix the problem of numerical instability which occurs when the sum of the stencil_coefficients is not
        exactly 0.
        https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb
        Stack Overflow answer on this issue.
        http://epubs.siam.org/doi/pdf/10.1137/S0036144596322507 - Modified Fornberg Algorithm
    =#
    _C = C[:, end]
    if order != 0
        _C[div(N, 2)+1] -= sum(_C)
    end
    return _C
end

index(i::Int, N::Int) = i + div(N, 2) + 1

struct DerivativeOperator{T<:Real,N,Wind,T2,S1,S2,S3,T3,F} <:
       AbstractDerivativeOperator{T}
    derivative_order::Int
    approximation_order::Int
    dx::T2
    len::Int
    stencil_length::Int
    stencil_coefs::S1
    boundary_stencil_length::Int
    boundary_point_count::Int
    low_boundary_coefs::S2
    high_boundary_coefs::S3
    offside::Int
    coefficients::T3
    coeff_func::F
end

struct CenteredDifference{N} end

function CenteredDifference{N}(derivative_order::Int,
    approximation_order::Int, dx::T,
    len::Int, coeff_func=1) where {T<:Real,N}
    @assert approximation_order > 1 "approximation_order must be greater than 1."
    stencil_length = derivative_order + approximation_order - 1 +
                     (derivative_order + approximation_order) % 2
    boundary_stencil_length = derivative_order + approximation_order
    dummy_x = (-div(stencil_length, 2)):div(stencil_length, 2)
    left_boundary_x = 0:(boundary_stencil_length-1)
    right_boundary_x = reverse((-boundary_stencil_length+1):0)

    boundary_point_count = div(stencil_length, 2) - 1 # -1 due to the ghost point
    # Because it's a N x (N+2) operator, the last stencil on the sides are the [b,0,x,x,x,x] stencils, not the [0,x,x,x,x,x] stencils, since we're never solving for the derivative at the boundary point.
    deriv_spots = (-div(stencil_length, 2)+1):-1  # unused
    L_boundary_deriv_spots = left_boundary_x[2:div(stencil_length, 2)]
    R_boundary_deriv_spots = right_boundary_x[2:div(stencil_length, 2)]

    stencil_coefs = convert(SVector{stencil_length,T},
        (1 / dx^derivative_order) *
        calculate_weights(derivative_order, zero(T), dummy_x))
    _low_boundary_coefs = SVector{boundary_stencil_length,T}[convert(SVector{
            boundary_stencil_length,
            T},
        (1 /
         dx^derivative_order) *
        calculate_weights(derivative_order,
            oneunit(T) *
            x0,
            left_boundary_x))
                                                             for x0 in L_boundary_deriv_spots]
    low_boundary_coefs = convert(SVector{boundary_point_count}, _low_boundary_coefs)

    # _high_boundary_coefs    = SVector{boundary_stencil_length, T}[convert(SVector{boundary_stencil_length, T}, (1/dx^derivative_order) * calculate_weights(derivative_order, oneunit(T)*x0, reverse(right_boundary_x))) for x0 in R_boundary_deriv_spots]
    high_boundary_coefs = convert(SVector{boundary_point_count},
        reverse(map(reverse,
            _low_boundary_coefs * (-1)^derivative_order)))

    offside = 0

    coefficients = fill!(Vector{T}(undef, len), 0)

    compute_coeffs!(coeff_func, coefficients)

    DerivativeOperator{T,N,false,T,typeof(stencil_coefs),
        typeof(low_boundary_coefs),typeof(high_boundary_coefs),
        typeof(coefficients),
        typeof(coeff_func)}(derivative_order, approximation_order, dx, len,
        stencil_length,
        stencil_coefs,
        boundary_stencil_length,
        boundary_point_count,
        low_boundary_coefs,
        high_boundary_coefs, offside, coefficients,
        coeff_func)
end

function generate_coordinates(i::Int, stencil_x, dummy_x,
    dx::AbstractVector{T}) where {T<:Real}
    len = length(stencil_x)
    stencil_x .= stencil_x .* zero(T)
    for idx in 1:div(len, 2)
        shifted_idx1 = index(idx, len)
        shifted_idx2 = index(-idx, len)
        stencil_x[shifted_idx1] = stencil_x[shifted_idx1-1] + dx[i+idx-1]
        stencil_x[shifted_idx2] = stencil_x[shifted_idx2+1] - dx[i-idx]
    end
    return stencil_x
end

function CenteredDifference{N}(derivative_order::Int,
    approximation_order::Int, dx::AbstractVector{T},
    len::Int, coeff_func=1) where {T<:Real,N}
    stencil_length = derivative_order + approximation_order - 1 +
                     (derivative_order + approximation_order) % 2
    boundary_stencil_length = derivative_order + approximation_order
    stencil_x = zeros(T, stencil_length)
    boundary_point_count = div(stencil_length, 2) - 1# -1 due to the ghost point

    interior_x = (boundary_point_count+2):(len+1-boundary_point_count)
    dummy_x = (-div(stencil_length, 2)):(div(stencil_length, 2)-1)
    low_boundary_x = [zero(T); cumsum(dx[1:(boundary_stencil_length-1)])]
    high_boundary_x = cumsum(dx[(end-boundary_stencil_length+1):end])
    # Because it's a N x (N+2) operator, the last stencil on the sides are the [b,0,x,x,x,x] stencils, not the [0,x,x,x,x,x] stencils, since we're never solving for the derivative at the boundary point.
    deriv_spots = (-div(stencil_length, 2)+1):-1

    stencil_coefs = [convert(SVector{stencil_length,T},
        calculate_weights(derivative_order, zero(T),
            generate_coordinates(i, stencil_x, dummy_x,
                dx)))
                     for i in interior_x]
    _low_boundary_coefs = SVector{boundary_stencil_length,T}[convert(SVector{
            boundary_stencil_length,
            T},
        calculate_weights(derivative_order,
            low_boundary_x[i+1],
            low_boundary_x))
                                                             for i in 1:boundary_point_count]
    low_boundary_coefs = convert(SVector{boundary_point_count}, _low_boundary_coefs)
    _high_boundary_coefs = SVector{boundary_stencil_length,T}[convert(SVector{
            boundary_stencil_length,
            T},
        calculate_weights(derivative_order,
            high_boundary_x[end-i],
            high_boundary_x))
                                                              for i in boundary_point_count:-1:1]
    high_boundary_coefs = convert(SVector{boundary_point_count}, _high_boundary_coefs)

    offside = 0

    coefficients = zeros(T, len)

    compute_coeffs!(coeff_func, coefficients)

    DerivativeOperator{T,N,false,typeof(dx),typeof(stencil_coefs),
        typeof(low_boundary_coefs),typeof(high_boundary_coefs),
        typeof(coefficients),
        typeof(coeff_func)}(derivative_order, approximation_order, dx,
        len, stencil_length,
        stencil_coefs,
        boundary_stencil_length,
        boundary_point_count,
        low_boundary_coefs,
        high_boundary_coefs, offside, coefficients,
        coeff_func)
end

CenteredDifference(args...) = CenteredDifference{1}(args...)
use_winding(A::DerivativeOperator{T,N,Wind}) where {T,N,Wind} = Wind
diff_axis(A::DerivativeOperator{T,N}) where {T,N} = N
function ==(A1::DerivativeOperator, A2::DerivativeOperator)
    return all([eval(:($A1.$name == $A2.$name)) for name in fieldnames(DerivativeOperator)])
end

function compute_coeffs!(coeff_func::Number,
    current_coeffs::AbstractVector{T}) where {T<:Number}
    return current_coeffs .+= coeff_func
end

function compute_coeffs!(coeff_func::AbstractVector{T},
    current_coeffs::AbstractVector{T}) where {T<:Number}
    return current_coeffs[:] += coeff_func
end

function compute_coeffs!(coeff_func::Function,
    current_coeffs::AbstractVector{T}) where {T<:Number}
    if hasmethod(coeff_func, (Vector{T},))
        current_coeffs[:] = coeff_func(current_coeffs)
    else
        map!(coeff_func, current_coeffs, current_coeffs)
    end
    return current_coeffs
end

function compute_coeffs(coeff_func::Number,
    current_coeffs::AbstractVector{T}) where {T<:Number}
    current_coeffs .+ coeff_func
end
function compute_coeffs(coeff_func::AbstractVector{T},
    current_coeffs::AbstractVector{T}) where {T<:Number}
    coeff_func + current_coeffs
end

function compute_coeffs(coeff_func::Function,
    current_coeffs::AbstractVector{T}) where {T<:Number}
    if hasmethod(coeff_func, (Vector{T},))
        return coeff_func(current_coeffs)
    else
        return map(coeff_func, current_coeffs)
    end
end

function Base.copyto!(L::AbstractMatrix{T}, A::DerivativeOperator{T}, N::Int) where {T}
    bl = A.boundary_point_count
    stencil_length = A.stencil_length
    stencil_pivot = use_winding(A) ? (1 + stencil_length % 2) : div(stencil_length, 2)
    bstl = A.boundary_stencil_length

    coeff = A.coefficients
    get_coeff = if coeff isa AbstractVector
        i -> coeff[i]
    elseif coeff isa Number
        i -> coeff
    else
        i -> true
    end

    for i in 1:bl
        cur_coeff = get_coeff(i)
        cur_stencil = use_winding(A) && cur_coeff < 0 ? reverse(A.low_boundary_coefs[i]) :
                      A.low_boundary_coefs[i]
        L[i, 1:bstl] = cur_coeff * cur_stencil
    end

    for i in (bl+1):(N-bl)
        cur_coeff = get_coeff(i)
        stencil = eltype(A.stencil_coefs) <: AbstractVector ? A.stencil_coefs[i-bl] :
                  A.stencil_coefs
        cur_stencil = use_winding(A) && cur_coeff < 0 ? reverse(stencil) : stencil
        L[i, (i+1-stencil_pivot):(i-stencil_pivot+stencil_length)] = cur_coeff *
                                                                     cur_stencil
    end

    for i in (N-bl+1):N
        cur_coeff = get_coeff(i)
        cur_stencil = use_winding(A) && cur_coeff < 0 ?
                      reverse(A.high_boundary_coefs[i-N+bl]) :
                      A.high_boundary_coefs[i-N+bl]
        L[i, (N-bstl+3):(N+2)] = cur_coeff * cur_stencil
    end

    L
end

function LinearAlgebra.Array(A::DerivativeOperator{T}, N::Int=A.len) where {T}
    copyto!(zeros(T, N, N + 2), A, N)
end

function SparseArrays.SparseMatrixCSC(A::DerivativeOperator{T}, N::Int=A.len) where {T}
    bl = A.boundary_point_count
    stencil_length = A.stencil_length
    stencil_pivot = use_winding(A) ? (1 + stencil_length % 2) : div(stencil_length, 2)
    bstl = A.boundary_stencil_length

    coeff = A.coefficients
    get_coeff = if coeff isa AbstractVector
        i -> coeff[i]
    elseif coeff isa Number
        i -> coeff
    else
        i -> true
    end

    Is = Int[]
    Js = Int[]
    Vs = T[]

    nvalues = 2 * bl * bstl + (N - 2 * bl) * stencil_length
    sizehint!(Is, nvalues)
    sizehint!(Js, nvalues)
    sizehint!(Vs, nvalues)

    for i in 1:bl
        cur_coeff = get_coeff(i)
        cur_stencil = use_winding(A) && cur_coeff < 0 ? reverse(A.low_boundary_coefs[i]) :
                      A.low_boundary_coefs[i]
        append!(Is, ((i for j in 1:bstl)...))
        append!(Js, 1:bstl)
        append!(Vs, cur_coeff * cur_stencil)
    end

    for i in (bl+1):(N-bl)
        cur_coeff = get_coeff(i)
        stencil = eltype(A.stencil_coefs) <: AbstractVector ? A.stencil_coefs[i-bl] :
                  A.stencil_coefs
        cur_stencil = use_winding(A) && cur_coeff < 0 ? reverse(stencil) : stencil
        append!(Is, ((i for j in 1:stencil_length)...))
        append!(Js, (i+1-stencil_pivot):(i-stencil_pivot+stencil_length))
        append!(Vs, cur_coeff * cur_stencil)
    end

    for i in (N-bl+1):N
        cur_coeff = get_coeff(i)
        cur_stencil = use_winding(A) && cur_coeff < 0 ?
                      reverse(A.high_boundary_coefs[i-N+bl]) :
                      A.high_boundary_coefs[i-N+bl]
        append!(Is, ((i for j in (N-bstl+3):(N+2))...))
        append!(Js, (N-bstl+3):(N+2))
        append!(Vs, cur_coeff * cur_stencil)
    end

    # ensure efficient allocation
    @assert length(Is) == nvalues
    @assert length(Js) == nvalues
    @assert length(Vs) == nvalues

    return sparse(Is, Js, Vs, N, N + 2)
end

function SparseArrays.sparse(A::DerivativeOperator{T}, N::Int=A.len) where {T}
    SparseMatrixCSC(A, N)
end

function Base.copyto!(L::AbstractSparseArray{T}, A::DerivativeOperator{T}, N::Int) where {T}
    copyto!(L, sparse(A))
end

function BandedMatrices.BandedMatrix(A::DerivativeOperator{T}, N::Int=A.len) where {T}
    stencil_length = A.stencil_length
    bstl = A.boundary_stencil_length
    L = BandedMatrix{T}(Zeros(N, N + 2), (bstl - 3, bstl - 1))
    copyto!(L, A, N)
end

function Base.convert(::Type{Mat},
    A::DerivativeOperator) where {
    Mat<:Union{Array,SparseMatrixCSC,
        BandedMatrix}}
    Mat(A)
end

Base.convert(::Type{AbstractMatrix}, A::DerivativeOperator) = BandedMatrix(A)

export CenteredDifference

end
