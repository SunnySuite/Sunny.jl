function fourier_exchange_matrix!(Jq::Matrix{ComplexF64}, sys::System; q)
    @assert sys.mode in (:dipole, :dipole_uncorrected) "SU(N) mode not supported"
    q_supercell = to_supercell_rlu(sys, q)

    Na = nsites(sys)
    Jq .= 0
    Jq = reshape(Jq, 3, Na, 3, Na)

    for i in 1:Na
        for coupling in interactions_at(sys, i).pair
            (; isculled, bond, bilin, biquad) = coupling
            isculled && break
            iszero(biquad) || error("Biquadratic interactions not supported")

            (j, n) = bond_target(sys, i, bond)
            J = exp(2π * im * dot(q_supercell, n)) * Mat3(bilin*I)
            view(Jq, :, i, :, j) .+= J
            view(Jq, :, j, :, i) .+= J'
        end
    end

    if !isnothing(sys.ewald)
        Aq = ewald_interaction_tensor(sys, q_supercell)
        for i in 1:Na, j in 1:Na
            view(Jq, :, i, :, j) .+= ewald_coupling(Aq, sys, q_supercell, i, j)
        end
    end

    for i in 1:Na
        (; c2, c4, c6) = interactions_at(sys, i).onsite
        iszero(c4) && iszero(c6) || error("Single-ion anisotropy beyond quadratic order not supported")
        anisotropy = SA[c2[1]-c2[3]        c2[5] 0.5c2[2];
                              c2[5] -c2[1]-c2[3] 0.5c2[4];
                           0.5c2[2]     0.5c2[4]   2c2[3]]
        view(Jq, :, i, :, i) .+= 2 * anisotropy
    end

    Jq = reshape(Jq, 3Na, 3Na)
    @assert maxdiff(Jq, Jq') < 1e-12
    return hermitianpart!(Jq)
end

function fourier_exchange_matrix_sensitivity!(Jq, sys, label; q)
    @assert sys.mode in (:dipole, :dipole_uncorrected) "SU(N) mode not supported"
    @assert sys.dims == (1, 1, 1) "System must have only a single cell"
    q_supercell = to_supercell_rlu(sys, q)

    Na = natoms(sys.crystal)
    Jq .= 0
    Jq = reshape(Jq, 3, Na, 3, Na)

    # Ignore param.val because it appears linearly and vanishes under derivative
    (; onsites, pairs) = lookup_param(sys, label)

    for coupling in pairs
        (; isculled, bond, bilin, biquad) = coupling
        isculled && break
        iszero(biquad) || error("Biquadratic interactions not supported")

        (; i, j, n) = bond

        J = exp(2π * im * dot(q_supercell, n)) * Mat3(bilin*I)
        view(Jq, :, i, :, j) .+= J
        view(Jq, :, j, :, i) .+= J'
    end

    for (i, onsite_coupling) in onsites
        (; c2, c4, c6) = onsite_coupling
        iszero(c4) && iszero(c6) || error("Single-ion anisotropy beyond quadratic order not supported")
        anisotropy = SA[c2[1]-c2[3]        c2[5] 0.5c2[2];
                              c2[5] -c2[1]-c2[3] 0.5c2[4];
                           0.5c2[2]     0.5c2[4]   2c2[3]]
        view(Jq, :, i, :, i) .+= 2 * anisotropy
    end

    Jq = reshape(Jq, 3Na, 3Na)
    @assert maxdiff(Jq, Jq') < 1e-12
    return hermitianpart!(Jq)
end

function fourier_exchange_matrix(sys::System; q)
    Na = nsites(sys)
    Jq = zeros(ComplexF64, 3Na, 3Na)
    return fourier_exchange_matrix!(Jq, sys; q)
end

function fourier_exchange_matrix_sensitivity(sys::System, label; q)
    Na = natoms(sys.crystal)
    Jq = zeros(ComplexF64, 3Na, 3Na)
    return fourier_exchange_matrix_sensitivity!(Jq, sys, label; q)
end
