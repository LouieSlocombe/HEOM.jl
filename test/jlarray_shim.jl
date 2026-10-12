# JLArrays is the reference CPU "device" of the GPUArrays ecosystem, but it has no FFT.
# These plans apply FFTW to host copies, so a `JLArray` runs the device code paths of
# HEOM (no scalar indexing, no host/device mixing) without a GPU.
using JLArrays: JLArray
using GPUArraysCore: allowscalar
import FFTW
import LinearAlgebra

struct HostPlan{P}
    plan::P
end

FFTW.plan_rfft(x::JLArray, region; kwargs...) =
    HostPlan(FFTW.plan_rfft(Array(x), region; kwargs...))
FFTW.plan_brfft(x::JLArray, d::Integer, region; kwargs...) =
    HostPlan(FFTW.plan_brfft(Array(x), d, region; kwargs...))
function LinearAlgebra.mul!(y::JLArray, p::HostPlan, x::JLArray)
    copyto!(y, p.plan * Array(x))
    return y
end

allowscalar(false)
