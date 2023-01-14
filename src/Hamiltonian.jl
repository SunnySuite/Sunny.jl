# Hamiltonian model parameters and energy/force calculations

mutable struct HamiltonianCPU # -> Model
    ext_field       :: Vector{Vec3}
    anisos          :: Vector{SingleIonAnisotropy}

    # pairints        :: Vector{PairInteractions}
    heisenbergs     :: Vector{HeisenbergCPU}
    gen_coups       :: Vector{GeneralCouplingCPU}
    biq_coups       :: Vector{BiquadraticCPU}

    ewald           :: Union{Nothing, EwaldCPU}
end


function HamiltonianCPU(ints::Vector{<:AbstractInteraction}, crystal::Crystal, N)
    nb = nbasis(crystal)
    ext_field = zeros(Vec3, nb)
    anisos = fill(SingleIonAnisotropy(N), nb)

    heisenbergs = Vector{HeisenbergCPU}()
    gen_coups   = Vector{GeneralCouplingCPU}()
    biq_coups   = Vector{BiquadraticCPU}()

    ewald       = nothing

    for int in ints
        if isa(int, QuadraticInteraction)
            validate_quadratic_interaction(int, crystal)
            int_impl = convert_quadratic(int, crystal)
            if isa(int_impl, HeisenbergCPU)
                push!(heisenbergs, int_impl)
            elseif isa(int_impl, GeneralCouplingCPU)
                push!(gen_coups, int_impl)
            else
                error("Quadratic interaction failed to convert to known backend type.")
            end
        elseif isa(int, BiQuadraticInteraction)
            if N != 0
                println("FIXME: BiQuadratic interactions are INCORRECT in SU(N) mode.")
            end
            int_imp2 = convert_biquadratic(int, crystal)
            push!(biq_coups, int_imp2)
        else
            error("$(int) failed to convert to known backend type.")
        end
    end

    return HamiltonianCPU(ext_field, anisos, heisenbergs, gen_coups, biq_coups, ewald)
end


function validate_quadratic_interaction(int::QuadraticInteraction, crystal::Crystal)
    # Validate all interactions
    int_str = repr("text/plain", int)
    b = int.bond

    # Verify that both basis sites indexed actually exist
    if !(1 <= b.i <= nbasis(crystal)) || !(1 <= b.j <= nbasis(crystal))
        error("Provided interaction $int_str indexes a non-existent basis site.")
    end

    # Verify that the interactions are symmetry-consistent
    if !is_coupling_valid(crystal, b, int.J)
        println("Symmetry-violating interaction: $int_str.")
        println("Use `print_bond(crystal, $b)` for more information.")
        error("Interaction violates symmetry.")
    end

    # We previously checked whether any interactions wrapped the entire system.
    # This check is now disabled because it can be useful to set the system size
    # equal to the magnetic unit cell.
    #=
    bs = all_symmetry_related_bonds(crystal, b)
    for b′ in bs
        coeffs = crystal.lat_vecs \ displacement(crystal, b′)
        wrapping = [i for i = 1:3 if abs(coeffs[i]) >= latsize[i]/2 - 1e-10]
        if !isempty(wrapping)
            println("Warning: Interaction $int_str wraps the system along dimension(s) $wrapping.")
        end
    end
    =#
end


function propagate_anisotropies!(hamiltonian::HamiltonianCPU, cryst::Crystal, b::Int, op::DP.AbstractPolynomialLike, N::Int)
    iszero(op) && return 

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
    nb = size(dipoles, 4)

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

    for heisen in ℋ.heisenbergs
        E += energy(dipoles, heisen)
    end
    for gen_coup in ℋ.gen_coups
        E += energy(dipoles, gen_coup)
    end
    for biq_coup in ℋ.biq_coups
        E += energy(dipoles, biq_coup)
    end
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

    for heisen in ℋ.heisenbergs
        J = first(heisen.bondtable.data)
        for (bond, _) in sublat_bonds(heisen.bondtable, b)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += J * (s⋅s - s₀⋅s₀)
            else
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                ΔE += J * (Δs ⋅ sⱼ)
            end
        end
    end
    for gen_coup in ℋ.gen_coups
        for (bond, J) in sublat_bonds(gen_coup.bondtable, b)
            if bond.i == bond.j && iszero(bond.n)
                ΔE += dot(s, J, s) - dot(s₀, J, s₀)
            else
                sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
                ΔE += dot(Δs, J, sⱼ)
            end
        end
    end
    for biq_coup in ℋ.biq_coups
        for (bond, effB) in sublat_bonds(biq_coup.bondtable, b)
            # On-site biquadratic does not make sense
            @assert !(bond.i == bond.j && iszero(bond.n))

            sⱼ = dipoles[offsetc(cell, bond.n, latsize), bond.j]
            ΔE += effB * ((s ⋅ sⱼ)^2 - (s₀ ⋅ sⱼ)^2)
        end
    end

    if !isnothing(ℋ.ewald)
        ΔE += energy_delta(dipoles, ℋ.ewald, idx, s)
    end

    return ΔE
end



"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ ``.
"""
function set_forces!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ℋ::HamiltonianCPU)
    # KBTODO remove this hack!
    N = size(ℋ.anisos[1].matrep, 1)

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
    
    for heisen in ℋ.heisenbergs
        accum_force!(B, dipoles, heisen)
    end
    for gen_coup in ℋ.gen_coups
        accum_force!(B, dipoles, gen_coup)
    end
    for biq_coup in ℋ.biq_coups
        accum_force!(B, dipoles, biq_coup)
    end
    if !isnothing(ℋ.ewald)
        accum_force!(B, dipoles, ℋ.ewald)
    end
end








#=
struct Model
    crystal         :: Crystal
    units           :: PhysicalConsts            # Physical constants that determine unit system
    latsize         :: NTuple{3, Int}            # Size of lattice in unit cells

    Ns              :: Vector{Int}               # Dimension of local Hilbert space per basis site in unit cell
    gs              :: Vector{Mat3}              # g-tensor per basis site in the crystal unit cell

    # Interactions
    ext_field       :: Vector{Vec3}
    anisos          :: Vector{SingleIonAnisotropy}
    pairints        :: Vector{PairInteractions}
    ewald           :: Union{Nothing, EwaldCPU}
end

struct SpinConfig{N}
    rng              :: Random.Xoshiro

    dipoles          :: Array{Vec3, 4}            # Expected dipoles
    coherents        :: Array{CVec{N}, 4}         # Coherent states
    κs               :: Vector{Float64}           # Meaning depends on context:
                                                  #  N > 0 => Effective ket rescaling, Z → √κ Z
                                                  #  N = 0 => Dipole magnitude, |s| = κ

    # Buffers for dynamical integration
    dipole_buffers   :: Vector{Array{Vec3, 4}}
    coherent_buffers :: Vector{Array{CVec{N}, 4}}
end
=#
