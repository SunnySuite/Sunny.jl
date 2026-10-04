# The Gaussian state about which the bosons are expanded, and Wick's theorem
# about it. Every static correction in this directory is one operation: a boson
# polynomial P, whether the Hamiltonian, an observable or a spin component, is
# re-expanded about a Gaussian state. Split each operator as O_a = w_a + δO_a, a
# c-number displacement plus a fluctuation whose pairs contract to the connected
# correlations of the state. Wick's theorem then writes P as a sum over the ways
# of displacing some of its slots and contracting others in pairs, the remaining
# slots standing in a product normal ordered with respect to the state.
#
# THE ORDER RULE. A word of K bosons in the expansion of a spin polynomial of
# degree d carries s^(d-K/2), and a displacement is of order s^(-1/2), so a term
# with d displacements and c contractions left standing with k bosons is smaller
# than the k-boson word itself by s^-(d+c). Each displacement and each
# contraction is a factor of 1/s. Keeping at most one factor is therefore
# exactly the O(1/s) correction to every word: the mean field (one contraction
# of H₄), the tadpole (one of H₃), the quadratic Hamiltonian of the displaced
# structure (one displacement of H₃), the corrected one-magnon amplitude, and
# the depletion and tilt of the ordered moment. Keeping more would inject a
# partial set of the next order, which a protected zero mode amplifies to the
# order kept; see the truncation rule in ExpansionDipole.jl.

# ---- The vacuum ----

# Quadratic Hamiltonians of a vacuum, as the congruence factors described at
# `TabulatedVacuum(swt, Hs)` below
struct TabulatedVacuum
    As :: Array{ComplexF64, 5}
end

"""
    MagnonVacuum(swt::SpinWaveTheory[, correction])

The Gaussian state about which the 1/s corrections expand the bosons: the vacuum
of the LSWT quadratic Hamiltonian, modified by `correction`. Its
quasi-particles are the internal lines of every loop, and the mean fields are
contractions in it. Without a correction this is the ``1/s`` expansion proper.

A `correction` is either a list of quadratic boson monomials to add to the LSWT
Hamiltonian, such as the mean-field corrections return, or a [`TabulatedVacuum`](@ref),
as [`dressed_vacuum`](@ref) builds. It resums some class of higher-order terms,
e.g. renormalized energies on the internal lines. Since the full Hamiltonian
does not depend on this choice, the difference from LSWT is subtracted again as
a counterterm in the static self-energy, so the bare propagator of the Dyson
equation remains that of LSWT and results differ from the ``1/s`` expansion only
at the order neglected.

!!! tip "Relation to self-consistent schemes"

    Renormalized energies on the internal lines, with harmonic vertices and
    coherence factors, is the scheme of [Veillette, James and Essler, PRB **72**,
    134429 (2005)](https://doi.org/10.1103/PhysRevB.72.134429), whose Dyson
    equation likewise keeps the bare propagator of LSWT. Iterating
    [`dressed_vacuum`](@ref) generalizes it to many bands. A vacuum that
    reproduces its own mean field is self-consistent Hartree-Fock. Either scheme
    resums an incomplete class of higher-order terms, so it is uncontrolled and
    may gap a Goldstone mode.
"""
struct MagnonVacuum
    swt        :: SpinWaveTheory
    correction :: Union{Vector{BosonMonomial{2}}, TabulatedVacuum}
end

MagnonVacuum(swt::SpinWaveTheory) = MagnonVacuum(swt, BosonMonomial{2}[])

# Overwrites the LSWT Hamiltonian `H` at `q_reshaped` with that of the vacuum
apply_correction!(H, terms::Vector{BosonMonomial{2}}, q_reshaped) = accum_quadratic!(H, terms, q_reshaped)

# Quadratic Hamiltonian of the vacuum at a wavevector in reshaped RLU
function vacuum_hamiltonian!(H, vac::MagnonVacuum, q_reshaped)
    dynamical_matrix!(H, vac.swt, q_reshaped)
    apply_correction!(H, vac.correction, q_reshaped)
end

# Adds the counterterm of the vacuum, the LSWT Hamiltonian minus that of the
# vacuum, to the Nambu matrix δH
function accum_counterterm!(δH, vac::MagnonVacuum, q_reshaped)
    H = zeros(ComplexF64, size(δH))
    dynamical_matrix!(H, vac.swt, q_reshaped)
    δH .+= H
    vacuum_hamiltonian!(H, vac, q_reshaped)
    δH .-= H
end

# ---- Tabulated vacuum ----

"""
    TabulatedVacuum(swt::SpinWaveTheory, Hs)

A vacuum correction given by its quadratic Hamiltonians `Hs[:, :, i, j, k]` in
the boson basis of `dynamical_matrix`, sampled at the wavevectors `(Tuple(c) .-
1/2) ./ dims` of the reshaped Brillouin zone, for each `c` in
`CartesianIndices(dims)`, and interpolated between. Pass to [`MagnonVacuum`](@ref).

No band labels are involved, so band crossings are harmless. What is
interpolated is the positive definite ``A`` that carries the LSWT Hamiltonian
``H`` into each sample, ``A H A = H̃``, and an interpolated ``A`` keeps ``A H
A`` positive semidefinite with the zero modes of ``H``. The vacuum is therefore
stable, and its quasi-particles are gapless wherever those of LSWT are.
"""
function TabulatedVacuum(swt::SpinWaveTheory, Hs::Array{<: Number, 5})
    As = similar(Hs, ComplexF64)
    dims = size(Hs)[3:5]
    for c in CartesianIndices(dims)
        H = dynamical_matrix(swt, Vec3((Tuple(c) .- 1/2) ./ dims))
        # The positive definite square root, which the regularization of
        # `swt` keeps invertible at a soft mode of LSWT
        R = sqrt(Hermitian(H))
        view(As, :, :, c) .= R \ sqrt(Hermitian(R * view(Hs, :, :, c) * R)) / R
    end
    return TabulatedVacuum(As)
end

function apply_correction!(H, tab::TabulatedVacuum, q_reshaped)
    A = interpolate_periodic(tab.As, q_reshaped)
    H .= A' * H * A
    hermitianpart!(H)
end

# Multilinear interpolation, periodic in the reshaped Brillouin zone, of the
# matrices `As[:, :, c]` sampled at the half-offset grid points (c - 1/2)/dims
function interpolate_periodic(As, q_reshaped)
    dims = size(As)[3:5]
    x = ntuple(d -> mod(q_reshaped[d] * dims[d] - 1/2, dims[d]), 3)
    i = floor.(Int, x)
    f = x .- i
    ret = zeros(ComplexF64, size(As, 1), size(As, 2))
    for corner in CartesianIndices((0:1, 0:1, 0:1))
        δ = Tuple(corner)
        w = prod(d -> δ[d] == 1 ? f[d] : 1 - f[d], 1:3)
        iszero(w) && continue
        c = ntuple(d -> mod(i[d] + δ[d], dims[d]) + 1, 3)
        ret .+= w .* view(As, :, :, c...)
    end
    return ret
end

# Quasi-particles of the vacuum, as `bogoliubov!` returns them in `ws`. A
# quadratic Hamiltonian that is not positive definite has no vacuum, which is an
# instability of the structure when the Hamiltonian is the harmonic one.
function vacuum_bogoliubov!(ws::BogoliubovWorkspace, H, vac::MagnonVacuum, q_reshaped)
    vacuum_hamiltonian!(H, vac, q_reshaped)
    try
        return bogoliubov!(ws, H)
    catch err
        err isa PosDefException || rethrow()
        rethrow(InstabilityError("Quadratic Hamiltonian of the vacuum not positive definite at reshaped wavevector \
                                  $(vec3_to_string(q_reshaped))."))
    end
end

# The wavevector loop shared by every frequency-dependent momentum integral of
# this module. Calls `f(𝐩, w, T1, T2, ε1, ε2)` once per wavevector 𝐩 of
# `grid`, where T1 = T(𝐩) and T2 = T(𝐪-𝐩) diagonalize the vacuum at the two
# magnon lines being paired, ε1, ε2 are the signed energies that `bogoliubov!`
# returns with them, and `w` is the multiplicity by which the contribution is to
# be scaled. All of the arrays are overwritten on each iteration, so `f` must
# consume them before returning.
function foreach_magnon_pair(f, vac::MagnonVacuum, q_reshaped, grid::LoopGrid)
    L = nbands(vac.swt)
    H = zeros(ComplexF64, 2L, 2L)
    # Both magnon lines are live at once, so each gets its own workspace
    (ws1, ws2) = (BogoliubovWorkspace(L), BogoliubovWorkspace(L))
    for (p, w) in zip(grid.ps, grid.wts)
        ε1 = vacuum_bogoliubov!(ws1, H, vac, p)
        ε2 = vacuum_bogoliubov!(ws2, H, vac, q_reshaped - p)
        f(p, w, ws1.T, ws2.T, ε1, ε2)
    end
end

# ---- Contractions ----

# Canonical label (a, a′, Δ) for the correlation ⟨δO_a(𝐫) δO_{a′}(𝐫+Δ)⟩,
# together with a flag indicating that the stored value is to be conjugated.
# Taking the adjoint of a correlation reverses and bars it,
#
#     conj⟨O_a(𝐫) O_{a′}(𝐫+Δ)⟩ = ⟨O_{ā′}(𝐫) O_{ā}(𝐫-Δ)⟩,
#
# and keeping only one representative of each such pair is what makes a
# mean-field Hamiltonian exactly Hermitian, rather than Hermitian only to within
# the accuracy of the momentum-space integration. The cell offset is labeled by
# integers, which are exactly comparable, unlike the floating point `Vec3` used
# elsewhere.
function correlation_key(a, a′, Δ::Vec3, L)
    Δ′ = round.(Int, Tuple(Δ))
    @assert Vec3(Δ′) ≈ Δ
    k = (a, a′, Δ′)
    kadj = (nambu_conj(a′, L), nambu_conj(a, L), .-Δ′)
    return kadj < k ? (kadj, true) : (k, false)
end

# Connected two-point correlations of the vacuum of a quadratic Hamiltonian, at
# the labels needed to contract the given monomials, and callable as `g(a, a′,
# Δ)` in either convention of `correlation_key`. A displacement leaves connected
# correlations unchanged, so only the quadratic Hamiltonian matters here. The
# `noise` is that of the quadrature, against which identities that hold only for
# exact integrals are asserted.
struct Contractions
    L      :: Int
    index  :: Dict{Tuple{Int, Int, NTuple{3, Int}}, Int}
    values :: Vector{ComplexF64}
    noise  :: Float64
end

function (g::Contractions)(a, a′, Δ)
    (k, conjugate) = correlation_key(a, a′, Δ, g.L)
    v = g.values[g.index[k]]
    return conjugate ? conj(v) : v
end

# The same labels with new values, e.g. a damped self-consistent update.
Contractions(g::Contractions, values) = Contractions(g.L, g.index, values, g.noise)

# Contractions of every pair of slots in `terms`, any iterable of monomials, in
# the vacuum. In momentum space,
#
#     ⟨x_𝐪[a] x_{-𝐪}[a′]⟩ = Σ_{n ≤ L} T_𝐪[a,n] conj(T_𝐪[ā′,n]),
#
# where the identity x_{-𝐪}[ā′] = x_𝐪[a′]† avoids a second diagonalization at
# -𝐪. That matters because `bogoliubov!` fixes the phase of each band
# independently, so only expressions built from a single T are gauge invariant.
# Averaging over the Brillouin zone with the phase exp(-2πi 𝐪⋅Δ) gives the
# real-space result.
#
# All requested labels share one quadrature, so that consumers of a common state
# pay for a single pass over the zone. On the loop grid of the cubic
# self-energy, a Ward identity relates these averages to that integrand point by
# point, which keeps each Goldstone mode exactly gapless, grid by grid.
function contractions(vac::MagnonVacuum, terms, quad::BZQuadrature)
    L = nbands(vac.swt)
    index = Dict{Tuple{Int, Int, NTuple{3, Int}}, Int}()
    for (; as, ns) in terms, p in eachindex(as), q in p+1:lastindex(as)
        (k, _) = correlation_key(as[p], as[q], ns[q] - ns[p], L)
        get!(index, k, length(index) + 1)
    end
    keys = first.(sort!(collect(index); by=last))
    noise = quadrature_noise(quad)
    isempty(keys) && return Contractions(L, index, ComplexF64[], noise)

    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    values = bz_average(quad) do q_reshaped
        vacuum_bogoliubov!(ws, H, vac, q_reshaped)
        U = view(ws.T, :, 1:L)
        return ComplexF64[cis(-2π * dot(q_reshaped, Vec3(Δ))) * dot(view(U, nambu_conj(a′, L), :), view(U, a, :))
                          for (a, a′, Δ) in keys]
    end
    return Contractions(L, index, values, noise)
end

# ---- Wick's theorem ----

# Every way to re-expand K slots leaving the k slots `S` standing, displacing
# the slots `D` and contracting the pairs `P`, with at most `nfactors`
# displacements and contractions in all. Each slot set is increasing, which
# keeps standing and contracted slots in their operator order.
function wick_patterns(K, k, nfactors)
    matchings(r) = isempty(r) ? [Tuple{Int, Int}[]] :
        [[(r[1], r[j]); m] for j in 2:length(r) for m in matchings(r[[2:j-1; j+1:end]])]
    subsets(r) = [r[findall(digits(Bool, m; base=2, pad=length(r)))] for m in 0:2^length(r)-1]

    ret = Tuple{Vector{Int}, Vector{Int}, Vector{Tuple{Int, Int}}}[]
    for S in subsets(collect(1:K))
        length(S) == k || continue
        for D in subsets(setdiff(1:K, S))
            R = setdiff(1:K, S, D)
            iseven(length(R)) && length(D) + length(R) ÷ 2 <= nfactors || continue
            append!(ret, (S, D, P) for P in matchings(R))
        end
    end
    return ret
end

# The coefficients of the normal-ordered products of k fluctuations into which
# each monomial re-expands about the Gaussian state of displacement `w`,
# Nambu-packed, and contractions `g`, keeping at most `nfactors` displacements
# and contractions per term; see the order rule above. A contracted pair (p, q),
# p < q, contributes ⟨δO_p δO_q⟩ in operator order. Passing `nothing` for `w`
# means no displacement, and for `g` no contraction.
#
# As an operator, :δO_r δO_s: differs from δO_r δO_s by a constant only, so a
# quadratic result may be used as a correction to the quadratic Hamiltonian of
# the fluctuations; for k = 0 the result is the expectation value in the state.
function wick_reduce(terms::Vector{BosonMonomial{K}}, g, w, ::Val{k}; nfactors=1) where {K, k}
    ret = BosonMonomial{k}[]
    for (S, D, P) in wick_patterns(K, k, nfactors)
        (isnothing(w) && !isempty(D) || isnothing(g) && !isempty(P)) && continue
        for (; c, as, ns) in terms
            for d in D
                c *= w[as[d]]
            end
            for (p, q) in P
                c *= g(as[p], as[q], ns[q] - ns[p])
            end
            iszero(c) || push!(ret, BosonMonomial(c, ntuple(i -> as[S[i]], Val{k}()), ntuple(i -> ns[S[i]], Val{k}())))
        end
    end
    return merge_monomials(ret)
end

# Expectation value of the sum of the monomials
wick_expectation(terms, g, w; nfactors=1) = sum(t -> t.c, wick_reduce(terms, g, w, Val{0}(); nfactors); init=0.0im)

# Coefficients ℓ[a] of a one-boson operator Σ_a ℓ[a] Σ_𝐫 O_a(𝐫), given as
# monomials. Cell offsets drop out, the operator being summed over all cells.
# Hermiticity of the operator, ℓ[ā] = conj(ℓ[a]), is asserted to within `noise`:
# it holds only for exact correlations, so the tolerance is set by the accuracy
# of their momentum integrals, whereas an error in a contraction would appear at
# O(1). The tolerance is absolute as well as relative, because symmetry can make
# ℓ vanish identically, leaving only integration noise.
function nambu_vector(terms::Vector{BosonMonomial{1}}, L, noise)
    ℓ = zeros(ComplexF64, 2L)
    for (; c, as) in terms
        ℓ[as[1]] += c
    end
    @assert norm(ℓ[L+1:2L] - conj(ℓ[1:L])) < noise * max(norm(ℓ), 1)
    return ℓ
end
