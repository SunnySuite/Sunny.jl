struct StevensExpansionP{T<:AbstractFloat}
    c0   :: SVector{1, T}
    c2   :: SVector{5, T}
    c4   :: SVector{9, T}
    c6   :: SVector{13, T}
    kmax :: Int
end
function StevensExpansionP{T}(src::Sunny.StevensExpansion) where T
    StevensExpansionP{T}(SVector{1,T}(src.c0), SVector{5,T}(src.c2),
                         SVector{9,T}(src.c4), SVector{13,T}(src.c6), src.kmax)
end
StevensExpansionP{Float64}(src::Sunny.StevensExpansion) =
    StevensExpansionP{Float64}(src.c0, src.c2, src.c4, src.c6, src.kmax)

struct PairCouplingDeviceKA{T<:AbstractFloat}
    isculled :: Bool
    bond     :: Sunny.Bond
    bilin    :: Union{T, SMatrix{3, 3, T, 9}}
    biquad   :: Union{T, SMatrix{5, 5, T, 25}}
end

function PairCouplingDeviceKA{T}(host::Sunny.PairCoupling) where T
    bilin  = host.bilin  isa Real ? T(host.bilin)            : SMatrix{3,3,T,9}(host.bilin)
    biquad = host.biquad isa Real ? T(host.biquad)           : SMatrix{5,5,T,25}(host.biquad)
    PairCouplingDeviceKA{T}(host.isculled, host.bond, bilin, biquad)
end
PairCouplingDeviceKA(host::Sunny.PairCoupling) = PairCouplingDeviceKA{Float64}(host)

function Adapt.adapt_structure(to, p::PairCouplingDeviceKA{T}) where T
    PairCouplingDeviceKA{T}(p.isculled,
                             Adapt.adapt_structure(to, p.bond),
                             Adapt.adapt_structure(to, p.bilin),
                             Adapt.adapt_structure(to, p.biquad))
end

struct InteractionsDeviceKA{TSExpansion}
    onsite :: TSExpansion
end

function Adapt.adapt_structure(to, inter::InteractionsDeviceKA)
    onsite = Adapt.adapt_structure(to, inter.onsite)
    InteractionsDeviceKA(onsite)
end

struct EwaldDeviceKA{T, TDemag, TArrA}
    μ0_μB² :: T
    demag  :: TDemag
    A      :: TArrA
end

function Adapt.adapt_structure(to, sys::EwaldDeviceKA)
    μ0_μB² = Adapt.adapt_structure(to, sys.μ0_μB²)
    demag  = Adapt.adapt_structure(to, sys.demag)
    A      = Adapt.adapt_structure(to, sys.A)
    EwaldDeviceKA(μ0_μB², demag, A)
end

@enum SystemMode begin
    dipole
    dipole_uncorrected
    SUN
end

function mode_to_enum(sys::Sunny.System{N}) where N
    if sys.mode == :SUN
        return SUN
    elseif sys.mode == :dipole
        return dipole
    else
        @assert sys.mode == :dipole_uncorrected
        return dipole_uncorrected
    end
end

struct CrystalDeviceKA{TLatvecs, TRecipvecs, Tpositions}
    latvecs   :: TLatvecs
    recipvecs :: TRecipvecs
    positions :: Tpositions
end

function Adapt.adapt_structure(to, data::CrystalDeviceKA)
    latvecs   = Adapt.adapt_structure(to, data.latvecs)
    recipvecs = Adapt.adapt_structure(to, data.recipvecs)
    positions = Adapt.adapt_structure(to, data.positions)
    CrystalDeviceKA(latvecs, recipvecs, positions)
end

@inline Sunny.natoms(cryst::CrystalDeviceKA) = length(cryst.positions)

struct SystemDeviceKA{TCrystal, TArrField, TArrInt, TPairs, TIndices, TArrGs, TEwald, TDipole}
    original_crystal   :: TCrystal
    mode               :: SystemMode
    crystal            :: TCrystal
    dims               :: NTuple{3, Int}
    extfield           :: TArrField
    interactions_union :: TArrInt
    pairs              :: TPairs
    indices            :: TIndices
    gs                 :: TArrGs
    ewald              :: TEwald
    dipoles            :: TDipole
end

function Adapt.adapt_structure(to, sys::SystemDeviceKA)
    original_crystal   = Adapt.adapt_structure(to, sys.original_crystal)
    crystal            = Adapt.adapt_structure(to, sys.crystal)
    dims               = Adapt.adapt_structure(to, sys.dims)
    extfield           = Adapt.adapt_structure(to, sys.extfield)
    interactions_union = Adapt.adapt_structure(to, sys.interactions_union)
    pairs              = Adapt.adapt_structure(to, sys.pairs)
    indices            = Adapt.adapt_structure(to, sys.indices)
    gs                 = Adapt.adapt_structure(to, sys.gs)
    ewald              = Adapt.adapt_structure(to, sys.ewald)
    dipoles            = Adapt.adapt_structure(to, sys.dipoles)
    SystemDeviceKA(original_crystal, sys.mode, crystal, dims, extfield,
                   interactions_union, pairs, indices, gs, ewald, dipoles)
end
