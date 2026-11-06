# The long-range dipole-dipole interaction, truncated to the quadratic boson
# Hamiltonian of linear spin wave theory, as `QuadraticTable` is for the other
# couplings. An ordered pair of magnetic moments μᵢ = -gᵢ 𝐒ᵢ contributes (μᵢ A
# μⱼ)/2 to the energy, and the reversed pair (j, i) supplies the symmetric half.
# For an entangled system the moments are those of the bare sites of each unit.
# The quadratic Hamiltonian has two parts. The quadratic words of each spin,
# weighted by the molecular field of all the others, are local and independent
# of 𝐪, and are compiled into the `QuadraticTable` with the other couplings,
# including its monomial list `swt.monomials.quadratic`. The linear words of
# pairs of spins couple every pair of sites, and are applied as
#
#     H_E(𝐪) = T(𝐪)† Ã(𝐪) T(𝐪).
#
# Here T maps a Nambu vector x to a complex 3-vector on each bare site i, the
# linear part of gᵢ 𝐒ᵢ. A component a of x refers to a boson of the unit that
# contains i, and one 3-vector per site serves both halves of the Nambu space.
# Between sites, Ã(𝐪) couples by the Ewald tensor Aᵢⱼ(𝐪) at the cell
# displacement Δ = 𝐜ⱼ - 𝐜ᵢ, wrapped into the supercell with the winding phase
# exp(-2πi 𝐪⋅w), w = fld(Δ, dims), as in `ewald_tensor_at`. For an entangled
# system, the sites are the bare sites of the uncontracted system, and 𝐜ᵢ is
# the unwrapped cell of the unit. This phase convention is a gauge choice.
#
# The dense matrix is built by `accum_ewald!`. `SpinWaveTheoryLarge` applies the
# same operator by FFT, through an `EwaldOperator`.
#
# A `QuadraticEwald` holds only this pair term, which is not among the monomials
# of `swt.monomials`. The cubic and quartic terms of the dipole-dipole
# interaction are not built at all, so the 1/s corrections refuse it.
struct QuadraticEwald
    cache  :: EwaldTensorCache
    # The entries ptr[i]:ptr[i+1]-1 of bare site i give the Nambu labels and
    # coefficients of T
    ptr    :: Vector{Int32}
    labels :: Vector{Int32}
    coefs  :: Vector{CVec3}
    # The atom and cell 𝐜ᵢ of each bare site
    atoms  :: Vector{Int}
    cells  :: Vector{NTuple{3, Int}}
end

# The Ewald pair term about the local `frames` of the classical structure, and
# the quadratic monomials of its molecular field, or `nothing` and no monomials
# without dipole-dipole interactions
function quadratic_ewald(sys::System, frames)
    usys = uncontracted_system(sys)
    isnothing(usys.ewald) && return (nothing, BosonMonomial{2}[])
    sites = [(i, bi, c) for i in 1:nsites(sys) for (bi, c) in bare_indices_and_cells(sys, i)]

    # The molecular field 𝐡ᵢ = ∂E/∂𝐒ᵢ of the Ewald energy E, at the classical
    # dipoles of `frames`, the words of no boson
    dipoles = zeros(Vec3, size(usys.dipoles))
    for (i, bi, _) in sites
        dipoles[bi] = Vec3(ntuple(3) do α
            real(sum(m -> m.c, spin_monomials(sys, frames, α, i, Val{0}(); bi); init=0.0im))
        end)
    end
    h = zeros(Vec3, size(usys.dipoles))
    accum_ewald_grad!(h, dipoles, usys)

    field = BosonMonomial{2}[]
    entries = [Dict{Int32, CVec3}() for _ in 1:nsites(usys)]
    atoms = vec([to_atom(site) for site in eachsite(usys)])
    cells = Vector{NTuple{3, Int}}(undef, nsites(usys))
    for (i, bi, c) in sites
        cells[bi] = c .- 1
        for α in 1:3
            for m in spin_monomials(sys, frames, α, i, Val{2}(); bi)
                push!(field, BosonMonomial(h[bi][α] * m.c, m.as, m.ns))
            end
            for m in spin_monomials(sys, frames, α, i, Val{1}(); bi)
                a = Int32(only(m.as))
                entries[bi][a] = get(entries[bi], a, zero(CVec3)) + usys.gs[bi][:, α] * m.c
            end
        end
    end
    ptr = Int32[1; 1 .+ cumsum(length.(entries))]
    labels = [a for d in entries for a in keys(d)]
    coefs = [w for d in entries for w in values(d)]

    (; cache) = usys.ewald
    E = QuadraticEwald(cache, ptr, labels, coefs, atoms, cells)
    return (E, field)
end

# H += T† Ã(𝐪) T, the Ewald pair term at the wavevector `q_supercell`
function accum_ewald!(H, E::QuadraticEwald, q_supercell)
    (; cache, ptr, labels, coefs, atoms, cells) = E
    q = Vec3(q_supercell)
    A = ewald_interaction_tensor(cache, q)
    for i in eachindex(cells), j in eachindex(cells)
        Aᵢⱼ = ewald_tensor_at(A, cache.dims, q, atoms[i], atoms[j], cells[j] .- cells[i])
        for e in ptr[i]:ptr[i+1]-1, f in ptr[j]:ptr[j+1]-1
            H[labels[e], labels[f]] += dot(coefs[e], Aᵢⱼ * coefs[f])
        end
    end
    return H
end

# A rough estimate of the operator norm of the pair term over all 𝐪: the
# dipole-dipole energy μ₀μ_B² s |g|² / r³ of the largest moments at the nearest
# distance r between two sites. It overestimates when one close pair dominates.
function ewald_scale(usys::System)
    (; μ0_μB², cryst) = usys.ewald.cache
    m² = maximum(i -> usys.κs[i] * opnorm(usys.gs[i])^2, eachindex(usys.κs))
    # Each atom has a periodic image at the shortest lattice vector
    max_dist = minimum(norm, eachcol(cryst.latvecs))
    r = minimum(cryst.positions) do pt
        (idxs, offsets) = all_offsets_within_distance(cryst.latvecs, cryst.positions, pt; min_dist=1e-8, max_dist)
        minimum(zip(idxs, offsets)) do (i, n)
            norm(cryst.latvecs * (cryst.positions[i] + n - pt))
        end
    end
    return μ0_μB² * m² / r^3
end
