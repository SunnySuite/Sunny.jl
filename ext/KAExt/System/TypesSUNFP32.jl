struct SystemDeviceSUNF32KA{
    TBondJ, TBondN, TIndices,
    TWOutSelf, TWOutCross, TWOutAnom,
    TWInSelf, TWInCross, TWInAnom, TWOnsite,
}
    Ns          :: Int
    Na          :: Int

    bond_j_out  :: TBondJ
    bond_n_out  :: TBondN
    indices_out :: TIndices

    bond_i_in   :: TBondJ
    bond_n_in   :: TBondN
    indices_in  :: TIndices

    W_out_self  :: TWOutSelf
    W_out_cross :: TWOutCross
    W_out_anom  :: TWOutAnom
    W_in_self   :: TWInSelf
    W_in_cross  :: TWInCross
    W_in_anom   :: TWInAnom
    W_onsite    :: TWOnsite
end

function SystemDeviceSUNF32KA(host::Sunny.System, backend)
    @assert host.mode == :SUN "SystemDeviceSUNF32KA requires :SUN mode"
    @assert isnothing(host.ewald) "Ewald interactions not supported for SU(N) GPU path"
    @assert host.interactions_union isa Vector "expected homogeneous Vector{Interactions}"
    @assert host.dims == (1, 1, 1) "expected SWT-reshaped sys with dims=(1,1,1); pass swt.sys"

    N  = host.Ns[1]
    Nf = N - 1
    Na = length(host.interactions_union)

    nout = zeros(Int, Na)
    nin  = zeros(Int, Na)
    for (me, int) in enumerate(host.interactions_union)
        for pair in int.pair
            pair.isculled && continue
            nout[me] += 1
            nin[pair.bond.j] += 1
        end
    end

    indices_out_h = Vector{Int32}(undef, Na + 1)
    indices_in_h  = Vector{Int32}(undef, Na + 1)
    indices_out_h[1] = Int32(1)
    indices_in_h[1]  = Int32(1)
    for me in 1:Na
        indices_out_h[me + 1] = Int32(indices_out_h[me] + nout[me])
        indices_in_h[me + 1]  = Int32(indices_in_h[me]  + nin[me])
    end
    total_out = Int(indices_out_h[Na + 1]) - 1
    total_in  = Int(indices_in_h[Na + 1])  - 1

    bond_j_out_h = Vector{Int32}(undef, total_out)
    bond_n_out_h = Vector{SVector{3, Int32}}(undef, total_out)
    bond_i_in_h  = Vector{Int32}(undef, total_in)
    bond_n_in_h  = Vector{SVector{3, Int32}}(undef, total_in)

    d_out = max(total_out, 1)
    d_in  = max(total_in,  1)
    W_out_self_h  = zeros(ComplexF64, Nf, Nf, d_out)
    W_out_cross_h = zeros(ComplexF64, Nf, Nf, d_out)
    W_out_anom_h  = zeros(ComplexF64, Nf, Nf, d_out)
    W_in_self_h   = zeros(ComplexF64, Nf, Nf, d_in)
    W_in_cross_h  = zeros(ComplexF64, Nf, Nf, d_in)
    W_in_anom_h   = zeros(ComplexF64, Nf, Nf, d_in)
    W_onsite_h    = zeros(ComplexF64, Nf, Nf, Na)

    write_out = Vector{Int32}(indices_out_h[1:Na])
    write_in  = Vector{Int32}(indices_in_h[1:Na])

    for (me, int) in enumerate(host.interactions_union)

        op  = int.onsite
        ONN = op[N, N]
        for n in 1:Nf, f in 1:Nf
            W_onsite_h[f, n, me] = op[f, n] - (f == n ? ONN : zero(ONN))
        end

        for pair in int.pair
            pair.isculled && continue
            j      = pair.bond.j
            n_bond = SVector{3, Int32}(pair.bond.n)

            oidx = Int(write_out[me])
            bond_j_out_h[oidx] = Int32(j)
            bond_n_out_h[oidx] = n_bond

            for (Ak, Bk) in pair.general.data
                ANN = Ak[N, N]
                BNN = Bk[N, N]
                for n in 1:Nf
                    BNn = Bk[N, n]
                    BnN = Bk[n, N]
                    for f in 1:Nf
                        AfN = Ak[f, N]
                        W_out_self_h[f, n, oidx]  += (Ak[f, n] - (f == n ? ANN : zero(ANN))) * BNN
                        W_out_cross_h[f, n, oidx] += AfN * BNn
                        W_out_anom_h[f, n, oidx]  += AfN * BnN
                    end
                end
            end
            write_out[me] += Int32(1)

            iidx = Int(write_in[j])
            bond_i_in_h[iidx] = Int32(me)
            bond_n_in_h[iidx] = n_bond

            for (Ak, Bk) in pair.general.data
                ANN = Ak[N, N]
                BNN = Bk[N, N]
                for f in 1:Nf
                    BfN = Bk[f, N]
                    for n in 1:Nf
                        W_in_self_h[f, n, iidx]  += ANN * (Bk[f, n] - (f == n ? BNN : zero(BNN)))
                        W_in_cross_h[f, n, iidx] += Ak[N, n] * BfN
                        W_in_anom_h[f, n, iidx]  += Ak[n, N] * BfN
                    end
                end
            end
            write_in[j] += Int32(1)
        end
    end

    @assert all(Int(write_out[me]) == Int(indices_out_h[me + 1]) for me in 1:Na) "Outgoing write_pos sanity"
    @assert all(Int(write_in[j])   == Int(indices_in_h[j + 1])   for j  in 1:Na) "Incoming write_pos sanity"

    W_out_self_d  = _backend_array(backend, ComplexF32.(W_out_self_h))
    W_out_cross_d = _backend_array(backend, ComplexF32.(W_out_cross_h))
    W_out_anom_d  = _backend_array(backend, ComplexF32.(W_out_anom_h))
    W_in_self_d   = _backend_array(backend, ComplexF32.(W_in_self_h))
    W_in_cross_d  = _backend_array(backend, ComplexF32.(W_in_cross_h))
    W_in_anom_d   = _backend_array(backend, ComplexF32.(W_in_anom_h))
    W_onsite_d    = _backend_array(backend, ComplexF32.(W_onsite_h))
    bond_j_out_d  = _backend_array(backend, bond_j_out_h)
    bond_n_out_d  = _backend_array(backend, bond_n_out_h)
    indices_out_d = _backend_array(backend, indices_out_h)
    bond_i_in_d   = _backend_array(backend, bond_i_in_h)
    bond_n_in_d   = _backend_array(backend, bond_n_in_h)
    indices_in_d  = _backend_array(backend, indices_in_h)

    return SystemDeviceSUNF32KA(N, Na,
        bond_j_out_d, bond_n_out_d, indices_out_d,
        bond_i_in_d, bond_n_in_d, indices_in_d,
        W_out_self_d, W_out_cross_d, W_out_anom_d,
        W_in_self_d, W_in_cross_d, W_in_anom_d,
        W_onsite_d)
end

function Adapt.adapt_structure(to, sys::SystemDeviceSUNF32KA)
    bond_j_out  = Adapt.adapt_structure(to, sys.bond_j_out)
    bond_n_out  = Adapt.adapt_structure(to, sys.bond_n_out)
    indices_out = Adapt.adapt_structure(to, sys.indices_out)
    bond_i_in   = Adapt.adapt_structure(to, sys.bond_i_in)
    bond_n_in   = Adapt.adapt_structure(to, sys.bond_n_in)
    indices_in  = Adapt.adapt_structure(to, sys.indices_in)
    W_out_self  = Adapt.adapt_structure(to, sys.W_out_self)
    W_out_cross = Adapt.adapt_structure(to, sys.W_out_cross)
    W_out_anom  = Adapt.adapt_structure(to, sys.W_out_anom)
    W_in_self   = Adapt.adapt_structure(to, sys.W_in_self)
    W_in_cross  = Adapt.adapt_structure(to, sys.W_in_cross)
    W_in_anom   = Adapt.adapt_structure(to, sys.W_in_anom)
    W_onsite    = Adapt.adapt_structure(to, sys.W_onsite)
    SystemDeviceSUNF32KA(sys.Ns, sys.Na,
        bond_j_out, bond_n_out, indices_out,
        bond_i_in, bond_n_in, indices_in,
        W_out_self, W_out_cross, W_out_anom,
        W_in_self, W_in_cross, W_in_anom,
        W_onsite)
end

@inline Sunny.natoms(sys::SystemDeviceSUNF32KA) = sys.Na
