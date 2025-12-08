# Pair couplings are counted only once per bond
struct PairCouplingDevice
    isculled :: Bool
    bond     :: Bond
    bilin    :: Union{Float64, Sunny.Mat3} # Bilinear
    biquad   :: Union{Float64, Sunny.Mat5} # Biquadratic
end

PairCouplingDevice(host::Sunny.PairCoupling) = PairCouplingDevice(host.isculled, host.bond, host.bilin, host.biquad)

function Adapt.adapt_structure(to, sys::PairCouplingDevice)
    isculled = Adapt.adapt_structure(to, sys.isculled)
    bond = Adapt.adapt_structure(to, sys.bond)
    bilin = Adapt.adapt_structure(to, sys.bilin)
    biquad = Adapt.adapt_structure(to, sys.biquad)
    PairCouplingDevice(isculled, bond, bilin, biquad)
end

struct InteractionsDevice{TSExpansion, TPair}
    onsite  :: TSExpansion
    pair    :: TPair
end

InteractionsDevice(host::Sunny.Interactions, pair_idx) = InteractionsDevice(host.onsite, pair_idx)

function Adapt.adapt_structure(to, inter::InteractionsDevice)
    onsite = Adapt.adapt_structure(to, inter.onsite)
    pair = Adapt.adapt_structure(to, inter.pair)
    InteractionsDevice(onsite, pair)
end

struct EwaldDevice
    μ0_μB²   :: Float64               # Strength of dipole-dipole interactions
    demag    :: Sunny.Mat3                  # Demagnetization factor
    A        :: CUDA.CuArray{Sunny.Mat3, 5}        # Interaction matrices in real-space         [offset+1,i,j]
end

EwaldDevice(host::Sunny.Ewald) = EwaldDevice(μ0_μB², demag, CUDA.CuArray(host.A))
EwaldDevice(host::Nothing) = Nothing()

function Adapt.adapt_structure(to, sys::EwaldDevice)
    μ0_μB² = Adapt.adapt_structure(to, sys.μ0_μB²)
    demag = Adapt.adapt_structure(to, sys.demag)
    A = Adapt.adapt_structure(to, sys.A)
    EwaldDevice(μ0_μB², demag, A)
end

struct SystemDevice{TCrystal, TArrField, TArrInt, TPairs, TArrGs}
    crystal            :: TCrystal
    extfield           :: TArrField # External B field
    interactions_union :: TArrInt # Interactions
    pairs              :: TPairs
    gs                 :: TArrGs # g-tensor per atom in unit cell
    #ewald              :: Union{EwaldDevice, Nothing}
end

function SystemDevice(host::Sunny.System)
    crystal = CrystalDevice(host.crystal)
    extfield = CUDA.CuArray(host.extfield)
    gs = CUDA.CuArray(host.gs)
    pairs_h = PairCouplingDevice[]
    interactions_h = InteractionsDevice{Sunny.StevensExpansion, Pair{Int64, Int64}}[]

    for int in host.interactions_union
        first = length(pairs_h) + 1
        last = length(pairs_h) + length(int.pair)
        for pair in int.pair
            push!(pairs_h, PairCouplingDevice(pair))
        end
        a = InteractionsDevice(int, Pair(first, last))
        push!(interactions_h, a)
    end
    pairs_d = CuVector(pairs_h)
    interactions_d = CuVector(interactions_h)
    return SystemDevice(crystal, extfield, interactions_d, pairs_d, gs)
end

function Adapt.adapt_structure(to, sys::SystemDevice)
    crystal = Adapt.adapt_structure(to, sys.crystal)
    extfield = Adapt.adapt_structure(to, sys.extfield)
    interactions_union = Adapt.adapt_structure(to, sys.interactions_union)
    pairs = Adapt.adapt_structure(to, sys.pairs)
    gs = Adapt.adapt_structure(to, sys.gs)
    #ewald = Adapt.adapt_structure(to, sys.ewald)
    SystemDevice(crystal, extfield, interactions_union, pairs, gs)
end