struct IncomingBondCSRKA{TPairs, TIndices, TMap}
    incoming_pairs   :: TPairs
    incoming_indices :: TIndices

    pair_map         :: TMap
end

function Adapt.adapt_structure(to, csr::IncomingBondCSRKA)
    incoming_pairs   = Adapt.adapt_structure(to, csr.incoming_pairs)
    incoming_indices = Adapt.adapt_structure(to, csr.incoming_indices)
    pair_map         = Adapt.adapt_structure(to, csr.pair_map)
    IncomingBondCSRKA(incoming_pairs, incoming_indices, pair_map)
end

function build_empty_incoming_csr(backend)

    pairs_d    = _backend_array(backend, PairCouplingDeviceKA{Float32}[])
    indices_d  = _backend_array(backend, Int64[1])
    pair_map_d = _backend_array(backend, Int64[])
    return IncomingBondCSRKA(pairs_d, indices_d, pair_map_d)
end

function build_incoming_bond_csr_ka(host_sys::Sunny.System, backend, ::Type{T}=Float64) where T<:AbstractFloat
    @assert host_sys.mode in (:dipole, :dipole_uncorrected) "build_incoming_bond_csr_ka is dipole-only; use build_incoming_bond_csr_sun_ka for SU(N)"
    @assert host_sys.interactions_union isa Vector "expected SpinWaveTheory-reshaped homogeneous Vector{Interactions}"
    @assert host_sys.dims == (1, 1, 1) "expected SWT-reshaped host_sys with dims=(1,1,1); got dims=$(host_sys.dims). Pass swt.sys to build_incoming_bond_csr_ka, not the raw user system."

    Na = length(host_sys.interactions_union)

    counts = zeros(Int, Na)
    for int in host_sys.interactions_union
        for pair in int.pair
            pair.isculled && continue
            j = pair.bond.j
            @assert 1 <= j <= Na "Bond j index $(j) out of range [1, $Na]"
            counts[j] += 1
        end
    end

    incoming_indices_h = Vector{Int64}(undef, Na + 1)
    incoming_indices_h[1] = 1
    for j in 1:Na
        incoming_indices_h[j+1] = incoming_indices_h[j] + counts[j]
    end
    total_incoming = incoming_indices_h[Na+1] - 1

    incoming_pairs_h = Vector{PairCouplingDeviceKA{T}}(undef, total_incoming)

    write_pos = Vector{Int}(undef, Na)
    for j in 1:Na
        write_pos[j] = incoming_indices_h[j]
    end

    for int in host_sys.interactions_union
        for pair in int.pair
            pair.isculled && continue
            j = pair.bond.j
            incoming_pairs_h[write_pos[j]] = PairCouplingDeviceKA{T}(pair)
            write_pos[j] += 1
        end
    end

    @assert all(write_pos[j] == incoming_indices_h[j+1] for j in 1:Na) "IncomingBondCSRKA write_pos sanity check failed"

    pairs_d    = _backend_array(backend, incoming_pairs_h)
    indices_d  = _backend_array(backend, incoming_indices_h)
    pair_map_d = _backend_array(backend, Int64[])

    return IncomingBondCSRKA(pairs_d, indices_d, pair_map_d)
end

function build_incoming_bond_csr_sun_ka(host_sys::Sunny.System, backend)
    @assert host_sys.mode == :SUN "build_incoming_bond_csr_sun_ka requires :SUN mode"
    @assert host_sys.interactions_union isa Vector "expected SpinWaveTheory-reshaped homogeneous Vector{Interactions}"
    @assert host_sys.dims == (1, 1, 1) "expected SWT-reshaped host_sys with dims=(1,1,1); pass swt.sys not raw system"

    Na = length(host_sys.interactions_union)

    counts = zeros(Int, Na)
    for int in host_sys.interactions_union
        for pair in int.pair
            pair.isculled && continue
            j = pair.bond.j
            @assert 1 <= j <= Na "Bond j index $(j) out of range [1, $Na]"
            counts[j] += 1
        end
    end

    incoming_indices_h = Vector{Int64}(undef, Na + 1)
    incoming_indices_h[1] = 1
    for j in 1:Na
        incoming_indices_h[j+1] = incoming_indices_h[j] + counts[j]
    end
    total_incoming = incoming_indices_h[Na+1] - 1

    incoming_pairs_h = Vector{PairCouplingDeviceKA{Float64}}(undef, total_incoming)
    pair_map_h       = Vector{Int64}(undef, total_incoming)

    write_pos = Vector{Int}(undef, Na)
    for j in 1:Na
        write_pos[j] = incoming_indices_h[j]
    end

    outgoing_flat = 0
    for int in host_sys.interactions_union
        for pair in int.pair
            outgoing_flat += 1
            pair.isculled && continue
            j = pair.bond.j
            slot = write_pos[j]
            incoming_pairs_h[slot] = PairCouplingDeviceKA{Float64}(pair)
            pair_map_h[slot] = outgoing_flat
            write_pos[j] += 1
        end
    end

    @assert all(write_pos[j] == incoming_indices_h[j+1] for j in 1:Na) "IncomingBondCSRKA (SUN) write_pos sanity check failed"

    pairs_d    = _backend_array(backend, incoming_pairs_h)
    indices_d  = _backend_array(backend, incoming_indices_h)
    pair_map_d = _backend_array(backend, pair_map_h)

    return IncomingBondCSRKA(pairs_d, indices_d, pair_map_d)
end
