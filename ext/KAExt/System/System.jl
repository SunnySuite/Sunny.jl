function _backend_array(backend, host_arr::AbstractArray)
    T = eltype(host_arr)
    d = KernelAbstractions.allocate(backend, T, size(host_arr)...)
    copyto!(d, host_arr)
    return d
end

_precconv(x::Real,                  ::Type{T}) where T<:AbstractFloat = T(x)
_precconv(x::Complex,               ::Type{T}) where T<:AbstractFloat = Complex{T}(x)
_precconv(x::SArray,                ::Type{T}) where T<:AbstractFloat = map(xi -> _precconv(xi, T), x)
_precconv(x::Sunny.StevensExpansion, ::Type{T}) where T<:AbstractFloat = StevensExpansionP{T}(x)
_precconv(x::PairCouplingDeviceKA,  ::Type{T}) where T<:AbstractFloat =
    PairCouplingDeviceKA{T}(x.isculled, x.bond,
                             _precconv(x.bilin,  T),
                             _precconv(x.biquad, T))
_precconv(x::InteractionsDeviceKA,  ::Type{T}) where T<:AbstractFloat =
    InteractionsDeviceKA(_precconv(x.onsite, T))
_precconv(x,                        ::Type{T}) where T<:AbstractFloat = x

function _backend_array(backend, host_arr::AbstractArray, ::Type{T}) where T<:AbstractFloat
    T == Float64 && return _backend_array(backend, host_arr)
    converted = map(x -> _precconv(x, T), host_arr)
    ET = eltype(converted)
    d = KernelAbstractions.allocate(backend, ET, size(converted)...)
    copyto!(d, converted)
    return d
end

function CrystalDeviceKA(host::Sunny.Crystal, backend, ::Type{T}=Float64) where T<:AbstractFloat
    positions_d = _backend_array(backend, host.positions, T)
    return CrystalDeviceKA(_precconv(host.latvecs, T), _precconv(host.recipvecs, T), positions_d)
end

function EwaldDeviceKA(host::Sunny.Ewald, backend, ::Type{T}=Float64) where T<:AbstractFloat
    A_d = _backend_array(backend, host.A, T)
    return EwaldDeviceKA(T(host.μ0_μB²), _precconv(host.demag, T), A_d)
end

EwaldDeviceKA(::Nothing, backend, ::Type{T}=Float64) where T<:AbstractFloat = nothing

function SystemDeviceKA(host::Sunny.System, backend, ::Type{T}=Float64) where T<:AbstractFloat
    @assert host.mode in (:dipole, :dipole_uncorrected)
    original_crystal = CrystalDeviceKA(Sunny.orig_crystal(host), backend, T)
    mode_enum = mode_to_enum(host)
    @assert mode_enum in (dipole, dipole_uncorrected)
    crystal  = CrystalDeviceKA(host.crystal, backend, T)
    dims     = host.dims
    extfield = _backend_array(backend, host.extfield, T)
    ewald    = EwaldDeviceKA(host.ewald, backend, T)
    gs       = _backend_array(backend, host.gs, T)
    dipoles  = _backend_array(backend, host.dipoles, T)

    indices_h = Int64[]
    push!(indices_h, 1)
    pairs_h = PairCouplingDeviceKA{T}[]
    interactions_h = InteractionsDeviceKA{StevensExpansionP{T}}[]

    for int in host.interactions_union
        push!(indices_h, indices_h[end] + length(int.pair))
        for pair in int.pair
            push!(pairs_h, PairCouplingDeviceKA{T}(pair))
        end
        push!(interactions_h, InteractionsDeviceKA(StevensExpansionP{T}(int.onsite)))
    end
    pairs_d        = _backend_array(backend, pairs_h)
    indices_d      = _backend_array(backend, indices_h)
    interactions_d = _backend_array(backend, interactions_h)

    return SystemDeviceKA(original_crystal, mode_enum, crystal, dims, extfield,
                          interactions_d, pairs_d, indices_d, gs, ewald, dipoles)
end

function Sunny.global_position(sys::SystemDeviceKA, site)
    r = sys.crystal.positions[site[4]] + Sunny.Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

@inline eachsite(sys::SystemDeviceKA) = CartesianIndices(size(sys.dipoles))
Sunny.nsites(sys::SystemDeviceKA) = length(eachsite(sys))