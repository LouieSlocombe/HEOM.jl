"""
    wigneranimation(sol; indices = eachindex(sol.u), kwargs...)
    wigneranimation(states, grid; times = nothing, indices = eachindex(states), kwargs...)

Record signed Wigner heatmaps as a `Plots.Animation`. Load `Plots` first. Export with
`gif(animation, "wigner.gif"; fps = 20)` or `mp4(animation, "wigner.mp4"; fps = 20)`.

`states` is a nonempty vector of real, finite matrices on `grid`. `indices` selects
saved states in playback order; repeated and reversed indices are allowed. Only
selected states are validated and used to determine the symmetric colour limits,
which stay fixed throughout the animation (zero data uses `(-1, 1)`). Position and
momentum axes also stay fixed. Samples are neither clipped nor renormalised.

The solution form reads the grid and saved times from a Wigner–Moyal solution.
For explicit states, optional `times` must be a finite real vector with one entry
per state. Titles show `t = ...` when times are available, otherwise `Frame ...`.
Standard Plots attributes, including `clims`, `title` and `size`, override defaults.

Each selected state produces one frame; no time interpolation is performed. GIF/MP4
playback gives all frames equal duration, so solve with uniform `saveat` for playback
proportional to simulation time. Frame PNGs are stored in Plots' temporary directory.
"""
function wigneranimation end

"""
    marginalanimation(sol; indices = eachindex(sol.u), kwargs...)
    marginalanimation(states, grid; times = nothing, indices = eachindex(states), kwargs...)

Record position and momentum densities in two panels as a `Plots.Animation`.
Load `Plots` first, then export with `gif` or `mp4`. Densities use the same quadrature
as [`marginalplot`](@ref) and retain their original normalisation and negative values.

Each panel's density limits include zero and all selected frames, with 5% padding;
an identically zero panel uses `(-1, 1)`. This keeps scales fixed during playback.
Standard Plots attributes override defaults; for example, `ylims = (-0.1, 1.0)`
sets both panels' limits.

State selection, validation, time labels and frame timing follow
[`wigneranimation`](@ref). Use uniform saved times for constant-speed playback.
"""
function marginalanimation end
