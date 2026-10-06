# Stationary preparation for the small double-well example. Requires HEOM,
# LinearAlgebra, LinearSolve, OrdinaryDiffEqVerner and SciMLBase.successful_retcode
# in the including script.
# Dense single-ADO LU factors use about 1.6 GB at 49² points and 35 auxiliaries.
# This is an example-specific accelerator, not a general equilibrium solver.

struct TraceConstrainedHEOM{O} <: AbstractMatrix{Float64}
    operator::O
    seed::Vector{Float64}
    dimensions::NTuple{3,Int}
end

Base.size(A::TraceConstrainedHEOM) = (length(A.seed), length(A.seed))

function LinearAlgebra.mul!(out::AbstractVector, A::TraceConstrainedHEOM, x::AbstractVector)
    U, dU = reshape(x, A.dimensions), reshape(out, A.dimensions)
    heom!(dU, U, A.operator, 0.0)
    integral = phase_space_integral(physical_wigner(U), A.operator.grid)
    @. out = -out + integral * A.seed
    return out
end

struct DoubleWellBlockPreconditioner{F} <: AbstractMatrix{Float64}
    factors::Vector{F}
    block_size::Int
end

function Base.size(P::DoubleWellBlockPreconditioner)
    n = P.block_size * length(P.factors)
    return (n, n)
end

function LinearAlgebra.ldiv!(
    out::AbstractVector,
    P::DoubleWellBlockPreconditioner,
    x::AbstractVector,
)
    output, input = reshape(out, P.block_size, :), reshape(x, P.block_size, :)
    for a in eachindex(P.factors)
        ldiv!(@view(output[:, a]), P.factors[a], @view(input[:, a]))
    end
    return out
end

function double_well_equilibrium(seed, op; root_shift = 0.01)
    all(isodd, size(op.grid)) || throw(
        ArgumentError(
            "stationary preparation requires odd grid sizes to avoid Nyquist null modes",
        ),
    )
    op.hamiltonian isa HEOM.SpectralWignerMoyal ||
        throw(ArgumentError("this preparation helper requires spectral discretisation"))
    isfinite(root_shift) && root_shift > 0 ||
        throw(ArgumentError("root_shift must be finite and positive"))
    # This also validates a static operator, finite input and positive root trace.
    initial = equilibrate(seed, (0.0, 0.0), op, Vern7())
    U = initial.hierarchy ./ initial.restart.norm
    W = physical_wigner(U)
    nq, np = size(op.grid)
    points = nq * np
    threads = BLAS.get_num_threads()
    try
        BLAS.set_num_threads(1)
        # Assemble only the single-ADO diagonal block from the exact spectral
        # operator, including its counterterm and bath terminator diffusion.
        L = zeros(points, points)
        input, output = zeros(nq, np), zeros(nq, np)
        for j in 1:points
            fill!(input, 0)
            input[j] = 1
            HEOM.apply_wigner_moyal!(output, input, op.hamiltonian)
            if !iszero(op.bath.diffusion)
                HEOM.heom_derivative!(op.scratch, input, op.hamiltonian, op.momentum.second)
                output .+= op.bath.diffusion .* op.scratch
            end
            L[:, j] = vec(output)
        end
        cache = Dict{Float64,LU{Float64,Matrix{Float64},Vector{Int}}}()
        factors = map(eachindex(op.damping)) do a
            damping = op.damping[a]
            get!(cache, damping) do
                # The positive root shift regularises only the preconditioner.
                # Neither this shift nor the block approximation changes HEOM.
                block = -L + max(damping, root_shift) * I
                if a == 1
                    block .+= vec(W) * fill(op.grid.dq * op.grid.dp, 1, points)
                end
                lu(block)
            end
        end
        P = DoubleWellBlockPreconditioner(factors, points)
        # (-L_HEOM + |seed><trace|) U = seed enforces both stationarity and
        # unit root trace. The complete generator remains matrix-free.
        b = copy(vec(U))
        A = TraceConstrainedHEOM(op, b, size(U))
        linear = solve(
            LinearSolve.LinearProblem(A, b),
            KrylovJL_GMRES(; gmres_restart = 200, timemax = 120.0);
            Pr = P,
            maxiters = 5000,
            abstol = 1e-10,
            reltol = 1e-10,
        )
        successful_retcode(linear) ||
            error("stationary preparation failed: $(linear.retcode)")
        result = equilibrate(
            reshape(linear.u, size(U)),
            (0.0, 0.0),
            op,
            Vern7();
            stationarity_abstol = 1e-9,
            stationarity_reltol = 1e-8,
        )
        result.converged || error("full-hierarchy stationarity check failed")
        abs(result.restart.norm - 1) < 1e-8 || error("equilibrium root trace is not one")
        return result
    finally
        BLAS.set_num_threads(threads)
    end
end
