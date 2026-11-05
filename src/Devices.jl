abstract type AbstractDevice{FT} end

"""
    CPU(FT=Float64)

A device for calculations on the host, with floating point precision `FT`. Pass
it as the `device` argument of [`SpinWaveTheory`](@ref) or
[`SpinWaveTheoryLarge`](@ref).

!!! warning "Reduced precision"

    Precision `Float32` is experimental. It is faster and uses less memory, but
    round-off of about ``10^{-7}`` relative to the spectral bandwidth requires a
    larger `regularization`, about ``10^{-6}`` times the bandwidth, which must be
    set explicitly. This can significantly modify the intensity of soft modes,
    e.g., near a Goldstone mode. For a spectrum known to be gapped,
    `regularization=0` is also allowed.
"""
struct CPU{FT} <: AbstractDevice{FT} end

# Calculations on an accelerator through a KernelAbstractions `backend`, e.g.,
# `CUDABackend()`. Not yet supported by any calculator.
struct GPU{FT, B} <: AbstractDevice{FT}
    backend :: B
end

function check_float_type(FT)
    FT in (Float32, Float64) || error("Floating point type must be Float32 or Float64")
    return FT
end

CPU(FT::Type=Float64) = CPU{check_float_type(FT)}()
GPU(backend, FT::Type=Float32) = GPU{check_float_type(FT), typeof(backend)}(backend)

float_type(::AbstractDevice{FT}) where FT = FT

Base.show(io::IO, ::CPU{FT}) where FT = print(io, "CPU($FT)")
Base.show(io::IO, d::GPU{FT}) where FT = print(io, "GPU($(d.backend), $FT)")

# Errors unless calculator `name` supports `device`
function check_device_supported(device, name)
    device isa AbstractDevice || error("Expected a device, e.g., CPU(Float32)")
    device isa CPU || error("$device is not yet supported by $name")
end
