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

struct InteractionsDevice{TSExpansion}
    onsite  :: TSExpansion
end

function Adapt.adapt_structure(to, inter::InteractionsDevice)
    onsite = Adapt.adapt_structure(to, inter.onsite)
    InteractionsDevice(onsite)
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

@enum SystemMode begin
    dipole
    dipole_uncorrected
    SUN
end

function mode_to_enum(sys::System{N}) where N
    if sys.mode == :SUN
        return SUN
    elseif sys.mode == :dipole
        return dipole
    else @assert sys.mode == :dipole_uncorrected
        return dipole_uncorrected
    end
end

struct SystemDevice{TCrystal, TArrField, TArrInt, TPairs, TIndices, TArrGs, TDipole}
    original_crystal   :: TCrystal
    mode               :: SystemMode
    crystal            :: TCrystal
    dims               :: NTuple{3, Int}
    extfield           :: TArrField # External B field
    interactions_union :: TArrInt # Interactions
    pairs              :: TPairs
    indices            :: TIndices
    gs                 :: TArrGs # g-tensor per atom in unit cell
    #ewald              :: Union{EwaldDevice, Nothing}
    dipoles            :: TDipole # Expected dipoles
end

function SystemDevice(host::Sunny.System)
    @assert host.mode in (:dipole, :dipole_uncorrected)
    @assert isnothing(host.ewald)
    original_crystal = CrystalDevice(Sunny.orig_crystal(host))
    mode = mode_to_enum(host)
    @assert mode in (dipole, dipole_uncorrected)
    crystal = CrystalDevice(host.crystal)
    dims = host.dims
    extfield = CUDA.CuArray(host.extfield)
    gs = CUDA.CuArray(host.gs)
    dipoles = CUDA.CuArray(host.dipoles)

    indices_h = Int64[]
    push!(indices_h, 1)
    pairs_h = PairCouplingDevice[]
    interactions_h = InteractionsDevice{Sunny.StevensExpansion}[]

    for int in host.interactions_union
        push!(indices_h,indices_h[end]+length(int.pair))
        for pair in int.pair
            push!(pairs_h, PairCouplingDevice(pair))
        end
        a = InteractionsDevice(int.onsite)
        push!(interactions_h, a)
    end
    pairs_d = CuVector(pairs_h)
    indices_d = CuVector(indices_h)
    interactions_d = CuVector(interactions_h)
    return SystemDevice(original_crystal, mode, crystal, dims, extfield, interactions_d, pairs_d, indices_d, gs, dipoles)
end

function Adapt.adapt_structure(to, sys::SystemDevice)
    original_crystal = Adapt.adapt_structure(to, sys.original_crystal)
    crystal = Adapt.adapt_structure(to, sys.crystal)
    dims = Adapt.adapt_structure(to, sys.dims)
    extfield = Adapt.adapt_structure(to, sys.extfield)
    interactions_union = Adapt.adapt_structure(to, sys.interactions_union)
    pairs = Adapt.adapt_structure(to, sys.pairs)
    indices = Adapt.adapt_structure(to, sys.indices)
    gs = Adapt.adapt_structure(to, sys.gs)
    #ewald = Adapt.adapt_structure(to, sys.ewald)
    dipoles = Adapt.adapt_structure(to, sys.dipoles)
    SystemDevice(original_crystal, sys.mode, crystal, dims, extfield, interactions_union, pairs, indices, gs, dipoles)
end

struct SystemDeviceSUN{TCrystal, TArrField, TPairs, TIndices, TOnsite, TArrGs, TDipole, TGeneral}
    original_crystal   :: TCrystal
    mode               :: SystemMode
    crystal            :: TCrystal
    dims               :: NTuple{3, Int}
    extfield           :: TArrField # External B field
    pairs              :: TPairs
    indices            :: TIndices
    onsite             :: TOnsite
    gs                 :: TArrGs # g-tensor per atom in unit cell
    #ewald              :: Union{EwaldDevice, Nothing}
    dipoles            :: TDipole # Expected dipoles
    Ns                 :: Int
    general            :: TGeneral
end

function SystemDeviceSUN(host::Sunny.System)
    @assert host.mode == :SUN
    @assert isnothing(host.ewald)
    original_crystal = CrystalDevice(Sunny.orig_crystal(host))
    mode = mode_to_enum(host)
    @assert mode == SUN
    crystal = CrystalDevice(host.crystal)
    dims = host.dims
    extfield = CUDA.CuArray(host.extfield)
    gs = CUDA.CuArray(host.gs)
    dipoles = CUDA.CuArray(host.dipoles)
    Ns = host.Ns[1]

    pair_length = 0
    data_4_len = 0
    data_3_len = 2
    (data_1_len, data_2_len) = size(host.interactions_union[begin].pair[begin].general.data[begin][1])
    for (i, int) in enumerate(host.interactions_union)
        pair_length += length(int.pair)
        for pair in int.pair
            data_4_len = max(data_1_len, length(pair.general.data))
            @assert data_1_len == size(pair.general.data[begin][1], 1)
            @assert data_2_len == size(pair.general.data[begin][1], 2)
        end
    end

    general_h = zeros(ComplexF64, data_1_len, data_2_len, data_3_len, data_4_len, pair_length)
    pairs_h = PairCouplingDevice[]
    indices_h = Int64[]
    push!(indices_h, 1)
    N = size(host.interactions_union[begin].onsite, 1)
    onsite_h = Array{ComplexF64}(undef, N, N, length(host.interactions_union))
    for (i, int) in enumerate(host.interactions_union)
        first = length(pairs_h) + 1
        push!(indices_h, indices_h[end]+length(int.pair))
        for (j, pair) in enumerate(int.pair)
            push!(pairs_h, PairCouplingDevice(pair))
            for (k, general_data) in enumerate(pair.general.data)
                view(general_h,:,:,1,k,first - 1 + j) .= general_data[1]
                view(general_h,:,:,2,k,first - 1 + j) .= general_data[2]
            end
        end
        view(onsite_h, :, :, i) .= int.onsite
    end
    pairs_d = CuVector(pairs_h)
    indices_d = CuVector(indices_h)
    onsite_d = CuArray(onsite_h)
    general_d = CuArray(general_h)
    return SystemDeviceSUN(original_crystal, SUN, crystal, dims, extfield, pairs_d, indices_d, onsite_d, gs, dipoles, Ns, general_d)
end

function Adapt.adapt_structure(to, sys::SystemDeviceSUN)
    original_crystal = Adapt.adapt_structure(to, sys.original_crystal)
    crystal = Adapt.adapt_structure(to, sys.crystal)
    dims = Adapt.adapt_structure(to, sys.dims)
    extfield = Adapt.adapt_structure(to, sys.extfield)
    pairs = Adapt.adapt_structure(to, sys.pairs)
    indices = Adapt.adapt_structure(to, sys.indices)
    onsite = Adapt.adapt_structure(to, sys.onsite)
    gs = Adapt.adapt_structure(to, sys.gs)
    #ewald = Adapt.adapt_structure(to, sys.ewald)
    dipoles = Adapt.adapt_structure(to, sys.dipoles)
    general = Adapt.adapt_structure(to, sys.general)
    SystemDeviceSUN(original_crystal, sys.mode, crystal, dims, extfield, pairs, indices, onsite, gs, dipoles, sys.Ns, general)
end
