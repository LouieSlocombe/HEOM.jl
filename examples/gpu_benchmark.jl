# Time one HEOM right-hand-side evaluation on the CPU and on a GPU.
#
# The per-member reference loop, the batched kernels on the CPU, and the batched kernels on a
# CUDA device in Float64 and Float32 are compared on the same hierarchy. Run it from an
# environment with HEOM, KernelAbstractions and CUDA:
#
#     julia --project=. examples/gpu_benchmark.jl
#
# Set the environment variables `HEOM_N` (grid points per axis), `HEOM_PADE` (Pade poles)
# and `HEOM_DEPTH` to scale the hierarchy; its size is printed. Only the right-hand side is
# timed, which is the cost of every explicit time step.

using HEOM, KernelAbstractions, CUDA, Printf

n = parse(Int, get(ENV, "HEOM_N", "64"))             # grid points per axis
pade = parse(Int, get(ENV, "HEOM_PADE", "3"))        # Pade poles
depth = parse(Int, get(ENV, "HEOM_DEPTH", "8"))      # hierarchy depth
nevals = 5                                           # right-hand-side evaluations per timing

grid = PhaseSpaceGrid((-8.0, 8.0), n, (-8.0, 8.0), n)
V = harmonic_potential(; mass = 1.0, omega = 1.0)
bath = drude_lorentz_pade_bath(; reorganization = 0.3, cutoff = 1.0, kT = 0.5, pade)

members = Int(hierarchy_size(length(bath.rates), depth))
@printf(
    "grid %d×%d, %d bath modes, depth %d: %d members, %.1f MB per Float64 state\n",
    n,
    n,
    length(bath.rates),
    depth,
    members,
    8 * n^2 * members / 1e6
)

# `time_rhs` returns seconds per evaluation, after a warm-up call and with the device
# synchronised so that asynchronous kernel launches are counted.
function time_rhs(op, U, sync)
    dU = similar(U)
    heom!(dU, U, op, zero(eltype(U)))
    sync()
    elapsed = @elapsed begin
        for _ in 1:nevals
            heom!(dU, U, op, zero(eltype(U)))
        end
        sync()
    end
    return elapsed / nevals
end

W0 = on_grid((q, p) -> exp(-((q - 1)^2 + p^2) / 2) / 2π, grid)
host = heom_operator(grid; mass = 1.0, potential = V, bath, depth, scaled = true)
U = heom_problem(W0, (0.0, 1.0), host).u0
reference = time_rhs(host, U, () -> nothing)
@printf("%-34s %10.3f ms  (1.0×)\n", "CPU, per-member loop", 1e3 * reference)

cpu = heom_operator(
    grid;
    mass = 1.0,
    potential = V,
    bath,
    depth,
    scaled = true,
    batched = true,
)
seconds = time_rhs(cpu, U, () -> nothing)
@printf(
    "%-34s %10.3f ms  (%.1f×)\n",
    "CPU, batched kernels",
    1e3 * seconds,
    reference / seconds
)

if CUDA.functional()
    CUDA.allowscalar(false)
    println("device: ", CUDA.name(CUDA.device()))
    for T in (Float64, Float32)
        gpu = heom_operator(
            grid;
            mass = 1.0,
            potential = V,
            bath,
            depth,
            scaled = true,
            device = CuArray,
            float_type = T,
        )
        Ug = heom_problem(W0, (0.0, 1.0), gpu).u0
        seconds = time_rhs(gpu, Ug, () -> CUDA.synchronize())
        @printf(
            "%-34s %10.3f ms  (%.1f×)\n",
            "GPU, batched kernels, $T",
            1e3 * seconds,
            reference / seconds
        )
    end
else
    println("no functional CUDA device; skipped the GPU timings")
end
