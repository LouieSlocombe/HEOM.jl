abstract type AbstractTimeDependentPotential end

"""
    TimeDependentPotential(V; derivative = nothing)

A real potential `V(q,t)`. `derivative(q,t)` may supply its partial time derivative;
otherwise ForwardDiff differentiates only the scalar potential, never FFT buffers.
A plain two-argument `V(q,t)` is also accepted by the operator constructors. Wrap
functions accepting both one and two arguments explicitly to select driven dynamics.
Spatial dual-number support is required for a truncated Moyal series; time dual-number
support is required only when a solver requests a time derivative and none was supplied.
For discontinuous drives, use solver `tstops` at the discontinuities.
"""
struct TimeDependentPotential{V,D} <: AbstractTimeDependentPotential
    potential::V
    derivative::D
end
TimeDependentPotential(V; derivative = nothing) = TimeDependentPotential(V, derivative)
(V::TimeDependentPotential)(q, t) = V.potential(q, t)

"""
    DrivenPotential(V0, field, dipole; field_derivative = nothing)

Separable potential `V0(q) - field(t)*dipole(q)`. Operators precompute the spatial
Moyal symbols (or finite-difference matrices), so each RHS evaluation only updates
one scalar field coefficient. The minus sign gives a force `+field(t)` for `dipole(q)=q`.
`field_derivative(t)` optionally supplies `d(field)/dt`; otherwise ForwardDiff is
used when an implicit solver requests the partial time derivative.
"""
struct DrivenPotential{V,E,M,D} <: AbstractTimeDependentPotential
    potential::V
    field::E
    dipole::M
    field_derivative::D
end
DrivenPotential(V0, field, dipole; field_derivative = nothing) =
    DrivenPotential(V0, field, dipole, field_derivative)
(V::DrivenPotential)(q, t) = V.potential(q) - V.field(t) * V.dipole(q)

normalize_potential(V::AbstractTimeDependentPotential) = V
normalize_potential(V) =
    !applicable(V, 0.0) && applicable(V, 0.0, 0.0) ? TimeDependentPotential(V) : V
is_time_dependent(V) = normalize_potential(V) isa AbstractTimeDependentPotential
is_time_dependent(op::AbstractWignerMoyal) = op.drive !== nothing
potential_at(V, t) = potential_at(normalize_potential(V), t, Val(true))
potential_at(V, t, ::Val{true}) = V
potential_at(V::AbstractTimeDependentPotential, t, ::Val{true}) = q -> V(q, t)

# Preserve the separable form when adding the time-independent bath counterterm.
add_counterterm(V, coefficient) =
    add_counterterm(normalize_potential(V), coefficient, Val(true))
add_counterterm(V, coefficient, ::Val{true}) = q -> V(q) + coefficient * q^2
add_counterterm(V::TimeDependentPotential, coefficient, ::Val{true}) =
    TimeDependentPotential((q, t) -> V(q, t) + coefficient * q^2; derivative = V.derivative)
add_counterterm(V::DrivenPotential, coefficient, ::Val{true}) = DrivenPotential(
    q -> V.potential(q) + coefficient * q^2,
    V.field,
    V.dipole;
    field_derivative = V.field_derivative,
)

struct SeparableMoyalDrive{P,A}
    potential::P
    baseline::A
    coupling::A
end
struct GeneralMoyalDrive{P}
    potential::P
end

function with_drive(op::SpectralWignerMoyal, drive)
    return SpectralWignerMoyal(
        op.grid,
        op.mass,
        op.hbar,
        op.moyal_terms,
        op.velocity,
        op.kinetic_symbol,
        op.potential_symbol,
        op.W,
        op.Wq,
        op.Wp,
        op.Tq,
        op.Tp,
        op.forward_q,
        op.backward_q,
        op.forward_p,
        op.backward_p,
        drive,
    )
end
with_drive(op::FiniteDifferenceWignerMoyal, drive) = FiniteDifferenceWignerMoyal(
    op.grid,
    op.mass,
    op.hbar,
    op.moyal_terms,
    op.order,
    op.matrix,
    drive,
)

function build_operator(d::Spectral, grid, mass, hbar, V::DrivenPotential, terms)
    op = build_operator(d, grid, mass, hbar, V.potential, terms)
    coupling =
        moyal_symbol(V.dipole, grid.q, wavenumbers(length(grid.p), grid.dp), hbar, terms)
    coupling ./= length(grid.p)
    check_finite(coupling)
    return with_drive(op, SeparableMoyalDrive(V, copy(op.potential_symbol), coupling))
end
function build_operator(d::FiniteDifference, grid, mass, hbar, V::DrivenPotential, terms)
    op = build_operator(d, grid, mass, hbar, V.potential, terms)
    coupling = fd_potential_matrix(op, V.dipole)
    # Retain the union of the spatial bands even when the field currently vanishes.
    baseline = op.matrix + 0.0 * coupling
    copyto!(op.matrix, baseline)
    return with_drive(op, SeparableMoyalDrive(V, baseline, coupling))
end
function build_operator(d::Spectral, grid, mass, hbar, V::TimeDependentPotential, terms)
    return with_drive(
        build_operator(d, grid, mass, hbar, q -> zero(q), terms),
        GeneralMoyalDrive(V),
    )
end
function build_operator(
    d::FiniteDifference,
    grid,
    mass,
    hbar,
    V::TimeDependentPotential,
    terms,
)
    return with_drive(
        build_operator(d, grid, mass, hbar, q -> zero(q), terms),
        GeneralMoyalDrive(V),
    )
end

function checked_drive_value(value)
    value isa Real && isfinite(value) || throw(
        ArgumentError(
            "drive and its time derivative must return finite real values, got $value",
        ),
    )
    converted = Float64(value)
    isfinite(converted) ||
        throw(ArgumentError("drive coefficient must be finite in Float64"))
    return converted
end

function fd_kinetic_matrix(op::FiniteDifferenceWignerMoyal)
    nq, np = size(op.grid)
    return kron(
        spdiagm(-op.grid.p ./ op.mass),
        periodic_difference_matrix(1, op.order, nq, op.grid.dq),
    )
end
function fd_potential_matrix(op::FiniteDifferenceWignerMoyal, V)
    nq, np = size(op.grid)
    derivatives = moyal_derivatives(V, op.grid.q, op.moyal_terms)
    L = SparseArrays.spzeros(nq * np, nq * np)
    for s in 0:(op.moyal_terms-1)
        Dp = periodic_difference_matrix(2s + 1, op.order, np, op.grid.dp)
        L += kron(Dp, spdiagm(moyal_coefficient(s, op.hbar) .* derivatives[:, s+1]))
    end
    all(isfinite, L.nzval) ||
        throw(ArgumentError("driven finite-difference coefficients must be finite"))
    return L
end

update_generator!(op::AbstractWignerMoyal, t) = update_generator!(op, op.drive, t)
update_generator!(op, ::Nothing, t) = nothing
function update_generator!(op::SpectralWignerMoyal, drive::SeparableMoyalDrive, t)
    E = checked_drive_value(drive.potential.field(t))
    @. op.potential_symbol = drive.baseline - E * drive.coupling
    check_finite(op.potential_symbol)
    return nothing
end
function update_generator!(op::FiniteDifferenceWignerMoyal, drive::SeparableMoyalDrive, t)
    E = checked_drive_value(drive.potential.field(t))
    copyto!(op.matrix, drive.baseline - E * drive.coupling)
    all(isfinite, op.matrix.nzval) ||
        throw(ArgumentError("driven generator must be finite"))
    return nothing
end
function update_generator!(op::SpectralWignerMoyal, drive::GeneralMoyalDrive, t)
    op.potential_symbol .=
        moyal_symbol(
            potential_at(drive.potential, t),
            op.grid.q,
            wavenumbers(length(op.grid.p), op.grid.dp),
            op.hbar,
            op.moyal_terms,
        ) ./ length(op.grid.p)
    check_finite(op.potential_symbol)
    return nothing
end
function update_generator!(op::FiniteDifferenceWignerMoyal, drive::GeneralMoyalDrive, t)
    copyto!(
        op.matrix,
        fd_kinetic_matrix(op) + fd_potential_matrix(op, potential_at(drive.potential, t)),
    )
    return nothing
end

function potential_time_derivative(V::TimeDependentPotential, q, t)
    value =
        V.derivative === nothing ? ForwardDiff.derivative(s -> V.potential(q, s), t) :
        V.derivative(q, t)
    value isa Real && isfinite(value) ||
        throw(ArgumentError("potential time derivative must be finite and real"))
    return value
end
function field_time_derivative(V::DrivenPotential, t)
    value =
        V.field_derivative === nothing ? ForwardDiff.derivative(V.field, t) :
        V.field_derivative(t)
    return checked_drive_value(value)
end

# Apply just the potential commutator, useful for impulsive linear response.
function potential_action!(out, W, op::SpectralWignerMoyal, t = 0.0)
    update_generator!(op, t)
    return apply_potential_symbol!(out, W, op, op.potential_symbol)
end
function apply_potential_symbol!(out, W, op::SpectralWignerMoyal, symbol)
    check_size(out, op.grid)
    check_size(W, op.grid)
    op.W .= W
    mul!(op.Wp, op.forward_p, op.W)
    op.Wp .*= symbol
    mul!(out, op.backward_p, op.Wp)
    return nothing
end
function potential_action!(out, W, op::FiniteDifferenceWignerMoyal, t = 0.0)
    update_generator!(op, t)
    check_size(out, op.grid)
    check_size(W, op.grid)
    mul!(vec(out), op.matrix - fd_kinetic_matrix(op), vec(W))
    return nothing
end

wigner_tgrad!(out, W, op::AbstractWignerMoyal, t) = wigner_tgrad!(out, W, op, op.drive, t)
wigner_tgrad!(out, W, op, ::Nothing, t) = (fill!(out, 0); nothing)
function wigner_tgrad!(out, W, op::SpectralWignerMoyal, drive::SeparableMoyalDrive, t)
    apply_potential_symbol!(out, W, op, drive.coupling)
    out .*= -field_time_derivative(drive.potential, t)
    return nothing
end
function wigner_tgrad!(
    out,
    W,
    op::FiniteDifferenceWignerMoyal,
    drive::SeparableMoyalDrive,
    t,
)
    mul!(vec(out), drive.coupling, vec(W))
    out .*= -field_time_derivative(drive.potential, t)
    return nothing
end
function wigner_tgrad!(out, W, op::SpectralWignerMoyal, drive::GeneralMoyalDrive, t)
    symbol =
        moyal_symbol(
            q -> potential_time_derivative(drive.potential, q, t),
            op.grid.q,
            wavenumbers(length(op.grid.p), op.grid.dp),
            op.hbar,
            op.moyal_terms,
        ) ./ length(op.grid.p)
    check_finite(symbol)
    return apply_potential_symbol!(out, W, op, symbol)
end
function wigner_tgrad!(out, W, op::FiniteDifferenceWignerMoyal, drive::GeneralMoyalDrive, t)
    matrix = fd_potential_matrix(op, q -> potential_time_derivative(drive.potential, q, t))
    mul!(vec(out), matrix, vec(W))
    return nothing
end
