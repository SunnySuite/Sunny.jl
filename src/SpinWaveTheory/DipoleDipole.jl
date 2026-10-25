# The long-range dipole-dipole interaction, added to the dynamical matrix `H`.
# It is dense in q and so not compiled into the quadratic table. An ordered pair
# of magnetic moments μᵢ = -gᵢ Sᵢ contributes (μᵢ A μⱼ)/2 to the energy, and the
# reversed pair (j, i) supplies the symmetric half. For an entangled system the
# moments are those of the bare sites of each unit.
function accum_ewald!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_supercell)
    (; sys) = swt
    q_supercell = Vec3(q_supercell)
    Aq = ewald_interaction_tensor(sys, q_supercell)
    A0 = ewald_interaction_tensor(sys, zero(Vec3))
    if sys.mode == :SUN
        accum_ewald_SUN!(H, swt, Aq, A0, q_supercell)
    else
        accum_ewald_dipole!(H, swt, Aq, A0, q_supercell)
    end
end

# The Fourier-space coupling gᵢᵀ Aᵢⱼ gⱼ between the spins of site i and site j,
# and all of its periodic images, where A = ewald_interaction_tensor(sys,
# q_supercell). Sites (i, j) are by linear index. The cell offset Δ determines
# both Aᵢⱼ and the supercell winding number w = fld(Δ, dims).
# 
# Typically Δ will be inferred from (i, j), but entangled systems require an
# override. Here, Bloch phases follow the cells of the units rather than the
# wrapped bare sites. This phase convention is a gauge choice for entangled
# systems and must be applied consistently.
function ewald_coupling(A, @nospecialize(sys::System), q_supercell, i, j, Δ=nothing)
    (; dims, gs) = sys
    (ci, ai) = to_cell_atom(eachsite(sys)[i])
    (cj, aj) = to_cell_atom(eachsite(sys)[j])
    Δ = something(Δ, cj .- ci)
    cell = CartesianIndex(mod.(Δ, dims) .+ 1)
    w = Vec3(fld.(Δ, dims))
    Aᵢⱼ = cis(-2π * dot(q_supercell, w)) * A[cell, ai, aj]
    return gs[i]' * Aᵢⱼ * gs[j]
end

# Each site i carries one Holstein-Primakoff boson. In the local frame of the
# ordered moment, Sˣ + iSʸ ≈ √2sᵢ b and Sᶻ = sᵢ - b†b. The transverse components
# of J(q) give hopping and pairing between the sites, and Jᶻᶻ at q = 0, the
# molecular field of the ordered moments, shifts each site by -sⱼ J₀ᶻᶻ.
function accum_ewald_dipole!(H, swt, Aq, A0, q_supercell)
    (; sys, data) = swt
    (; local_rotations, sqrtS) = data
    L = nbands(swt)
    H11 = view(H, 1:L, 1:L)
    H12 = view(H, 1:L, L+1:2L)
    H21 = view(H, L+1:2L, 1:L)
    H22 = view(H, L+1:2L, L+1:2L)

    for i in 1:L, j in 1:L
        Jq = ewald_coupling(Aq, sys, q_supercell, i, j) / 2
        J0 = ewald_coupling(A0, sys, zero(Vec3), i, j) / 2
        Ri, Rj = local_rotations[i], local_rotations[j]
        J = sqrtS[i] * sqrtS[j] * Ri' * Jq * Rj
        J0zz = (Ri' * J0 * Rj)[3, 3]

        Q⁻ = (J[1, 1] + J[2, 2] - im*(J[1, 2] - J[2, 1])) / 2
        Q⁺ = (J[1, 1] + J[2, 2] + im*(J[1, 2] - J[2, 1])) / 2
        H11[i, j] += Q⁻
        H11[j, i] += conj(Q⁻)
        H22[i, j] += Q⁺
        H22[j, i] += conj(Q⁺)

        P⁻ = (J[1, 1] - J[2, 2] - im*(J[1, 2] + J[2, 1])) / 2
        P⁺ = (J[1, 1] - J[2, 2] + im*(J[1, 2] + J[2, 1])) / 2
        H21[i, j] += P⁻
        H21[j, i] += conj(P⁺)
        H12[i, j] += P⁺
        H12[j, i] += conj(P⁻)

        si, sj = sqrtS[i]^2, sqrtS[j]^2
        H11[i, i] -= sj * J0zz
        H22[i, i] -= sj * J0zz
        H11[j, j] -= si * J0zz
        H22[j, j] -= si * J0zz
    end
end

# Each unit i carries N-1 bosons bₘ. In the local frame, where the condensate is
# basis state N, an operator O expands to bilinear order as
#
#     O = O[N,N] + Σₘ (O[m,N] b†ₘ + O[N,m] bₘ) + Σₘₙ (O[m,n] - δₘₙ O[N,N]) b†ₘ bₙ.
#
# A coupling J A⊗B between units i and j is bilinear in two ways. A quadratic
# term on one unit times the mean on the other stays on one unit, and takes J
# at q = 0. Linear terms on both units couple them, and take J at q.
function accum_ewald_SUN!(H, swt, Aq, A0, q_supercell)
    (; sys, data) = swt
    (; spin_ops) = data::SWTDataSUN
    usys = uncontracted_system(sys)
    N = sys.Ns[1]
    Na = nsites(sys)
    L = (N-1) * Na
    blockdims = (N-1, Na, N-1, Na)
    H11 = reshape(view(H, 1:L, 1:L), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), blockdims)

    # The bosons live on units, so the phase of each coupling follows the cell
    # displacement of the units, not that of the wrapped bare sites
    for i in 1:Na, j in 1:Na, (bi, ci) in bare_indices_and_cells(sys, i), (bj, cj) in bare_indices_and_cells(sys, j)
        Jq = ewald_coupling(Aq, usys, q_supercell, bi, bj, cj .- ci) / 2
        J0 = ewald_coupling(A0, usys, zero(Vec3), bi, bj, cj .- ci) / 2
        for α in 1:3, β in 1:3
            A, B = spin_ops[α, bi], spin_ops[β, bj]
            jq, j0 = Jq[α, β], J0[α, β]
            for m in 1:N-1, n in 1:N-1
                # Quadratic on one unit, mean on the other
                ci = (A[m, n] - δ(m, n) * A[N, N]) * B[N, N]
                cj = A[N, N] * (B[m, n] - δ(m, n) * B[N, N])
                H11[m, i, n, i] += ci * j0
                H22[n, i, m, i] += ci * j0
                H11[m, j, n, j] += cj * j0
                H22[n, j, m, j] += cj * j0
                # Hopping i → j and j → i, and pairing
                hij = A[m, N] * B[N, n]
                hji = A[N, m] * B[n, N]
                p = A[m, N] * B[n, N]
                H11[m, i, n, j] += hij * jq
                H22[n, j, m, i] += hij * conj(jq)
                H11[n, j, m, i] += hji * conj(jq)
                H22[m, i, n, j] += hji * jq
                H12[m, i, n, j] += p * jq
                H12[n, j, m, i] += p * conj(jq)
                H21[n, j, m, i] += conj(p) * conj(jq)
                H21[m, i, n, j] += conj(p) * jq
            end
        end
    end
end
