# Hamiltonian model parameters and energy/force calculations


mutable struct HamiltonianCPU # -> Model
    ext_field :: Vector{Vec3}
    anisos    :: Vector{SingleIonAnisotropy}
    pairints  :: Vector{PairInteractions}
    ewald     :: Union{Nothing, EwaldCPU}
end

function HamiltonianCPU(crystal::Crystal, N)
    nb = nbasis(crystal)
    ext_field = zeros(Vec3, nb)
    anisos = fill(SingleIonAnisotropy(N), nb)
    pairints = [PairInteractions() for _ in 1:nb]
    ewald = nothing

    return HamiltonianCPU(ext_field, anisos, pairints, ewald)
end


function validate_atom_range(cryst::Crystal, a::Int)
    (1 <= a <= nbasis(cryst)) || error("Atom index $a is out of range.")
end

function propagate_exchange!(hamiltonian::HamiltonianCPU, cryst::Crystal, J::Mat3, J_biq::Float64, bond::Bond)
    # Verify that atom indices are in range
    validate_atom_range(cryst, bond.i)
    validate_atom_range(cryst, bond.j)

    # Verify that exchange is symmetry-consistent
    if !is_coupling_valid(cryst, bond, J)
        println("Symmetry-violating exchange: $J.")
        println("Use `print_bond(crystal, $bond)` for more information.")
        error("Interaction violates symmetry.")
    end

    # Biquadratic interactions not yet supported in SU(N) mode
    N = size(hamiltonian.anisos[1].matrep, 1)
    if !iszero(J_biq) && N > 0
        error("Biquadratic interactions not yet supported in SU(N) mode.")
    end

    # Print a warning if an interaction already exists for bond
    (; heisen, exchng, biquad) = hamiltonian.pairints[bond.i]
    if any(x -> x[2] == bond, vcat(heisen, exchng, biquad))
        println("Warning: Overriding exchange for bond $bond.")
    end

    isheisen = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)

    for (i, (; heisen, exchng, biquad)) in enumerate(hamiltonian.pairints)
        for (bond′, J′) in zip(all_symmetry_related_couplings_for_atom(cryst, i, bond, J)...)
            # Remove any existing interactions for bond′
            matches_bond(x) = x[2] == bond′
            filter!(!matches_bond, heisen)
            filter!(!matches_bond, exchng)
            filter!(!matches_bond, biquad)

            # The energy or force calculation only needs to see each bond once.
            # Use a trick below to select half the bonds to cull.
            bond_delta = (bond′.j - bond′.i, bond′.n...)
            @assert bond_delta != (0, 0, 0, 0)
            isculled = bond_delta > (0, 0, 0, 0)

            if isheisen
                @assert J ≈ J′
                push!(heisen, (isculled, bond′, J′[1,1]))
            else
                push!(exchng, (isculled, bond′, J′))
            end

            if !iszero(J_biq)
                push!(biquad, (isculled, bond′, J_biq))
            end
        end

        # Sort interactions so that non-culled bonds appear first
        sort!(heisen, by=first)
        sort!(exchng, by=first)
        sort!(biquad, by=first)
    end
end


function propagate_anisotropies!(hamiltonian::HamiltonianCPU, cryst::Crystal, b::Int, op::DP.AbstractPolynomialLike, N::Int)
    iszero(op) && return 

    validate_atom_range(cryst, b)

    if !iszero(hamiltonian.anisos[b].op)
        println("Warning: Overriding anisotropy for atom $b.")
    end

    if !is_anisotropy_valid(cryst, b, op)
        println("Symmetry-violating anisotropy: $op.")
        println("Use `print_site(crystal, $b)` for more information.")
        error("Invalid anisotropy.")
    end

    for (b′, op′) in zip(all_symmetry_related_anisotropies(cryst, b, op)...)
        matrep = operator_to_matrix(op′; N)

        S = (N-1)/2
        c = operator_to_classical_stevens_coefficients(op′, S)
        all(iszero.(c[[1,3,5]])) || error("Odd-ordered dipole anisotropies not supported.")
        c2 = SVector{ 5}(c[2])
        c4 = SVector{ 9}(c[4])
        c6 = SVector{13}(c[6])
        kmax = max(!iszero(c2)*2, !iszero(c4)*4, !iszero(c6)*6)
        clsrep = ClassicalStevensExpansion(kmax, c2, c4, c6)

        hamiltonian.anisos[b′] = SingleIonAnisotropy(op′, matrep, clsrep)
    end
end

function energy(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, ℋ::HamiltonianCPU, κs::Vector{Float64}) where N
    E = 0.0
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    @inbounds for idx in CartesianIndices(dipoles)
        E -= ℋ.ext_field[idx[4]] ⋅ dipoles[idx]
    end

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        for idx in CartesianIndices(coherents)
            E_idx, _ = energy_and_gradient_for_classical_anisotropy(dipoles[idx], ℋ.anisos[idx[4]].clsrep)
            E += E_idx
        end
    else
        for idx in CartesianIndices(coherents)
            Λ = ℋ.anisos[idx[4]].matrep
            κ = κs[idx[4]]
            Z = coherents[idx]
            E += κ * real(Z' * Λ * Z)
        end
    end

    for (; heisen, exchng, biquad) in ℋ.pairints
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += J * dot(sᵢ, sⱼ)
            end
        end
        # Quadratic exchange
        for (culled, bond, J) in exchng
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += dot(sᵢ, J, sⱼ)
            end
        end
        # Biquadratic exchange
        for (culled, bond, J) in biquad
            culled && break
            for cell in CartesianIndices(latsize)
                sᵢ = dipoles[cell, bond.i]
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                E += J * dot(sᵢ, sⱼ)^2
            end
        end
    end

    # Long-range dipole-dipole
    if !isnothing(ℋ.ewald)
        E += energy(dipoles, ℋ.ewald)
    end
    return E
end


# Computes the change in energy for an update to spin state
function energy_local_delta(dipoles::Array{Vec3, 4}, coherents::Array{CVec{N}, 4}, ℋ::HamiltonianCPU, κs::Vector{Float64}, idx, s::Vec3, Z::CVec{N}) where N
    s₀ = dipoles[idx]
    Z₀ = coherents[idx]
    Δs = s - s₀
    ΔE = 0.0

    cell, b = splitidx(idx)
    latsize = size(dipoles)[1:3]

    # Zeeman coupling to external field
    ΔE -= ℋ.ext_field[b] ⋅ Δs

    # Single-ion anisotropy, dipole or SUN mode
    if N == 0
        E_new, _ = energy_and_gradient_for_classical_anisotropy(s, ℋ.anisos[b].clsrep)
        E_old, _ = energy_and_gradient_for_classical_anisotropy(s₀, ℋ.anisos[b].clsrep)
        ΔE += E_new - E_old
    else
        Λ = ℋ.anisos[b].matrep
        ΔE += κs[b] * real(dot(Z, Λ, Z) - dot(Z₀, Λ, Z₀))
    end

    (; heisen, exchng, biquad) = ℋ.pairints[b]
    # Heisenberg exchange
    for (_, bond, J) in heisen
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * (Δs ⋅ sⱼ)
    end
    # Quadratic exchange
    for (_, bond, J) in exchng
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += dot(Δs, J, sⱼ)
    end
    # Biquadratic exchange
    for (_, bond, J) in biquad
        sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
        ΔE += J * ((s ⋅ sⱼ)^2 - (s₀ ⋅ sⱼ)^2)
    end

    if !isnothing(ℋ.ewald)
        ΔE += energy_delta(dipoles, ℋ.ewald, idx, s)
    end

    return ΔE
end



# Updates B in-place to hold negative energy gradient, -dE/ds, for each spin.
function set_forces!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU)
    # KBTODO remove this hack!
    N = size(ℋ.anisos[1].matrep, 1)
    latsize = size(B)[1:3]

    fill!(B, zero(Vec3))

    # Zeeman coupling
    @inbounds for idx in CartesianIndices(dipoles)
        B[idx] += ℋ.ext_field[idx[4]]
    end

    # Single-ion anisotropy only contributes in dipole mode. In SU(N) mode, the
    # anisotropy matrix will be incorporated directly into ℌ.
    if N == 0
        for idx in CartesianIndices(dipoles)
            _, gradE = energy_and_gradient_for_classical_anisotropy(dipoles[idx], ℋ.anisos[idx[4]].clsrep)
            B[idx] -= gradE
        end
    end
    
    for (; heisen, exchng, biquad) in ℋ.pairints
        # Heisenberg exchange
        for (culled, bond, J) in heisen
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= J  * sⱼ
                B[cellⱼ, bond.j] -= J' * sᵢ
            end
        end
        # Quadratic exchange
        for (culled, bond, J) in exchng
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= J  * sⱼ
                B[cellⱼ, bond.j] -= J' * sᵢ
            end
        end
        # Biquadratic exchange
        for (culled, bond, J) in biquad
            culled && break
            for cellᵢ in CartesianIndices(latsize)
                cellⱼ = offsetc(cellᵢ, bond.n, latsize)
                sᵢ = dipoles[cellᵢ, bond.i]
                sⱼ = dipoles[cellⱼ, bond.j]
                B[cellᵢ, bond.i] -= 2J  * sⱼ * (sᵢ⋅sⱼ)
                B[cellⱼ, bond.j] -= 2J' * sᵢ * (sᵢ⋅sⱼ)
            end
        end
    end

    if !isnothing(ℋ.ewald)
        accum_force!(B, dipoles, ℋ.ewald)
    end
end

