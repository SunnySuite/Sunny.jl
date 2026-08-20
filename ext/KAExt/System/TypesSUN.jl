struct SystemDeviceSUNKA{TCrystal, TArrField, TPairs, TIndices, TOnsite, TArrGs, TDipole, TGeneral}
    original_crystal :: TCrystal
    crystal          :: TCrystal
    mode             :: SystemMode
    dims             :: NTuple{3, Int}
    extfield         :: TArrField
    pairs            :: TPairs
    indices          :: TIndices
    onsite           :: TOnsite
    gs               :: TArrGs
    dipoles          :: TDipole
    Ns               :: Int
    general          :: TGeneral
end

function SystemDeviceSUNKA(host::Sunny.System, backend)
    @assert host.mode == :SUN "SystemDeviceSUNKA requires :SUN mode system"
    @assert isnothing(host.ewald) "Ewald interactions not supported for SU(N) GPU path"
    @assert host.interactions_union isa Vector "expected homogeneous Vector{Interactions}"
    @assert host.dims == (1, 1, 1) "expected SWT-reshaped sys with dims=(1,1,1); pass swt.sys not raw system"

    original_crystal = CrystalDeviceKA(Sunny.orig_crystal(host), backend)
    crystal  = CrystalDeviceKA(host.crystal, backend)
    dims     = host.dims
    extfield = _backend_array(backend, host.extfield)
    gs       = _backend_array(backend, host.gs)
    dipoles  = _backend_array(backend, host.dipoles)
    Ns       = host.Ns[1]
    N        = Ns
    Na       = length(host.interactions_union)

    total_pairs = 0
    max_general = 0
    for int in host.interactions_union
        total_pairs += length(int.pair)
        for pair in int.pair
            max_general = max(max_general, length(pair.general.data))
        end
    end

    general_h = zeros(ComplexF64, N, N, 2, max(max_general, 1), max(total_pairs, 1))
    onsite_h  = zeros(ComplexF64, N, N, Na)
    pairs_h   = PairCouplingDeviceKA{Float64}[]
    indices_h = Int64[1]

    for (i, int) in enumerate(host.interactions_union)
        push!(indices_h, indices_h[end] + length(int.pair))
        for (j, pair) in enumerate(int.pair)
            push!(pairs_h, PairCouplingDeviceKA{Float64}(pair))
            flat_idx = indices_h[i] + j - 1
            for (k, (Ak, Bk)) in enumerate(pair.general.data)
                view(general_h, :, :, 1, k, flat_idx) .= Ak
                view(general_h, :, :, 2, k, flat_idx) .= Bk
            end
        end
        view(onsite_h, :, :, i) .= int.onsite
    end

    pairs_d   = _backend_array(backend, pairs_h)
    indices_d = _backend_array(backend, indices_h)
    onsite_d  = _backend_array(backend, onsite_h)
    general_d = _backend_array(backend, general_h)

    return SystemDeviceSUNKA(original_crystal, crystal, SUN, dims,
                              extfield, pairs_d, indices_d, onsite_d,
                              gs, dipoles, Ns, general_d)
end

function Adapt.adapt_structure(to, sys::SystemDeviceSUNKA)
    original_crystal = Adapt.adapt_structure(to, sys.original_crystal)
    crystal          = Adapt.adapt_structure(to, sys.crystal)
    extfield         = Adapt.adapt_structure(to, sys.extfield)
    pairs            = Adapt.adapt_structure(to, sys.pairs)
    indices          = Adapt.adapt_structure(to, sys.indices)
    onsite           = Adapt.adapt_structure(to, sys.onsite)
    gs               = Adapt.adapt_structure(to, sys.gs)
    dipoles          = Adapt.adapt_structure(to, sys.dipoles)
    general          = Adapt.adapt_structure(to, sys.general)
    SystemDeviceSUNKA(original_crystal, crystal, sys.mode, sys.dims,
                       extfield, pairs, indices, onsite, gs, dipoles, sys.Ns, general)
end

@inline Sunny.natoms(sys::SystemDeviceSUNKA) = Sunny.natoms(sys.crystal)
