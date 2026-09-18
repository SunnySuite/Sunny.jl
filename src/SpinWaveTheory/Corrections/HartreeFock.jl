# Mean-field (Hartree-Fock) treatment of the four-boson term H₄, which shifts the
# magnon dispersion at O(1/s) relative to LSWT. Conventions are collected in
# Corrections.jl and the monomial representation in Vertices.jl.

# Nambu index of the adjoint operator, exchanging b_i and b†_i.
nambu_conj(a, L) = mod1(a + L, 2L)

# Accumulates quadratic monomials into the Nambu matrix H, in the convention
# H₂ = (1/2) x†_𝐪 H_𝐪 x_𝐪 + const of `swt_hamiltonian_dipole!`. Fourier
# transforming a monomial as in `vertex!` gives
#
#     c Σ_𝐫 O_{a₁}(𝐫+𝐧₁) O_{a₂}(𝐫+𝐧₂) = Σ_𝐪 c φ x_𝐪[a₁] x_𝐪[ā₂]†,
#
# where φ = exp(2πi 𝐪⋅(𝐧₁-𝐧₂)) and ā = mod1(a+L, 2L) labels the adjoint. Up to a
# constant from the commutator this is c φ x†[ā₂] x[a₁]. Since Σ_𝐪 x†_𝐪[α] x_𝐪[β]
# is invariant under (α, β, φ) → (β̄, ᾱ, conj(φ)), the weight can be split evenly
# between the two forms; that is what feeds a single b†b monomial into both the
# (1,1) and the (2,2) Nambu block. Applied to the monomials of H₂ itself, the rule
# reproduces every entry written by `swt_hamiltonian_dipole!`.
#
# The rule is covariant under taking adjoints, so a monomial list describing a
# Hermitian operator yields a Hermitian matrix. Two monomials can nonetheless
# describe the same operator while being written differently, in which case
# Hermiticity holds only to the accuracy with which their coefficients were
# computed; `hermitianpart!` projects out that residual.
function accum_quadratic!(H, terms::Vector{BosonMonomial{2}}, q_reshaped)
    L = div(size(H, 1), 2)
    for (; c, as, ns) in terms
        φ = cis(2π * dot(q_reshaped, ns[1] - ns[2]))
        H[nambu_conj(as[2], L), as[1]] += c * φ
        H[nambu_conj(as[1], L), as[2]] += c * conj(φ)
    end
    hermitianpart!(H)
end

# The three ways to split four slots into two pairs, each pair keeping the
# original order of its slots.
const QUARTIC_PAIRINGS = ((1, 2, 3, 4), (1, 3, 2, 4), (1, 4, 2, 3))

# Canonical label (a, a′, Δ) for the ground-state correlation
# ⟨O_a(𝐫) O_{a′}(𝐫+Δ)⟩, together with a flag indicating that the stored value is
# to be conjugated. Taking the adjoint of a correlation reverses and bars it,
#
#     conj⟨O_a(𝐫) O_{a′}(𝐫+Δ)⟩ = ⟨O_{ā′}(𝐫) O_{ā}(𝐫-Δ)⟩,
#
# and keeping only one representative of each such pair is what makes the
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

# Canonical labels of every correlation needed to decouple the given monomials.
function correlation_keys(L, termss...)
    ret = Tuple{Int, Int, NTuple{3, Int}}[]
    for terms in termss, (; as, ns) in terms, p in eachindex(as), q in p+1:lastindex(as)
        (k, _) = correlation_key(as[p], as[q], ns[q] - ns[p], L)
        k in ret || push!(ret, k)
    end
    return ret
end

# Ground-state correlations ⟨O_a(𝐫) O_{a′}(𝐫+Δ)⟩ for each requested key. In
# momentum space,
#
#     ⟨x_𝐪[a] x_{-𝐪}[a′]⟩ = Σ_{n ≤ L} T_𝐪[a,n] conj(T_𝐪[ā′,n]),
#
# where the identity x_{-𝐪}[ā′] = x_𝐪[a′]† avoids a second diagonalization at -𝐪.
# That matters because `bogoliubov!` fixes the phase of each band independently,
# so only expressions built from a single T are gauge invariant. Averaging over
# the Brillouin zone with the phase exp(-2πi 𝐪⋅Δ) gives the real-space result.
#
# The extra quadratic terms `terms2` allow the correlations to be evaluated in the
# already corrected ground state, as required for self-consistency.
function nambu_correlations(swt::SpinWaveTheory, ckeys, terms2; opts...)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L)

    (gs, _) = hcubature((0, 0, 0), (1, 1, 1); opts...) do q
        q_reshaped = Vec3(q)
        dynamical_matrix!(H, swt, q_reshaped)
        accum_quadratic!(H, terms2, q_reshaped)
        bogoliubov!(T, H)
        U = view(T, :, 1:L)
        return [cis(-2π * dot(q_reshaped, Vec3(Δ))) * dot(view(U, nambu_conj(a′, L), :), view(U, a, :))
                for (a, a′, Δ) in ckeys]
    end

    return gs
end

# Wraps the output of `nambu_correlations` as a function of the labels
# (a, a′, Δ) in either of the two conventions of `correlation_key`.
function correlation_lookup(ckeys, gs, L)
    G = Dict(zip(ckeys, gs))
    return function (a, a′, Δ)
        (k, conjugate) = correlation_key(a, a′, Δ, L)
        return conjugate ? conj(G[k]) : G[k]
    end
end

# Wick decoupling of a four-boson term into a quadratic form plus a constant. For
# each pairing of the slots,
#
#     O_p O_q O_r O_s → ⟨O_p O_q⟩ O_r O_s + ⟨O_r O_s⟩ O_p O_q - ⟨O_p O_q⟩⟨O_r O_s⟩.
#
# Subtracting the constant once, rather than twice, is what makes the expectation
# value of the mean-field form agree with ⟨H₄⟩ itself.
function hartree_fock_decoupling(terms::Vector{BosonMonomial{4}}, g)
    terms2 = BosonMonomial{2}[]
    δE = zero(ComplexF64)

    for (; c, as, ns) in terms, (p, q, r, s) in QUARTIC_PAIRINGS
        g1 = g(as[p], as[q], ns[q] - ns[p])
        g2 = g(as[r], as[s], ns[s] - ns[r])
        push!(terms2, BosonMonomial(c * g1, (as[r], as[s]), (ns[r], ns[s])))
        push!(terms2, BosonMonomial(c * g2, (as[p], as[q]), (ns[p], ns[q])))
        δE -= c * g1 * g2
    end

    return (merge_monomials(terms2), δE)
end

"""
    hartree_fock_correction(swt::SpinWaveTheory; maxiters=1, tol=1e-8, damping=0, opts...)

Decouples the four-boson term of the Holstein-Primakoff expansion into a
mean-field correction to the quadratic (LSWT) Hamiltonian. The correction is
smaller than the LSWT Hamiltonian by a factor of order ``1/s``. Returns
`(; terms2, δE, iters)`, where `terms2` can be passed to
[`corrected_dispersion`](@ref) and `δE` is a correction to the energy per site.

With `maxiters=1` the mean fields are those of the uncorrected LSWT ground state.
Larger values iterate to self-consistency, stopping when the mean fields move by
less than `tol`. A nonzero `damping` in `[0, 1)` mixes in the previous iterate,
which can stabilize the iteration.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration.
"""
function hartree_fock_correction(swt::SpinWaveTheory; maxiters=1, tol=1e-8, damping=0, opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")
    check_corrections_supported(swt)

    L = nbands(swt)
    terms4 = quartic_monomials(swt)
    ckeys = correlation_keys(L, terms4)
    terms2 = BosonMonomial{2}[]
    gs = zeros(ComplexF64, length(ckeys))
    δE = zero(ComplexF64)
    iters = 0

    for iter in 1:maxiters
        gs′ = nambu_correlations(swt, ckeys, terms2; opts...)
        iter > 1 && (gs′ = damping*gs + (1-damping)*gs′)
        converged = iter > 1 && norm(gs′ - gs) < tol
        gs = gs′
        (terms2, δE) = hartree_fock_decoupling(terms4, correlation_lookup(ckeys, gs, L))
        iters = iter
        converged && break
    end

    @assert abs(imag(δE)) < 1e-9 * max(abs(δE), 1)
    return (; terms2, δE = real(δE) / nsites(uncontracted_system(swt.sys)), iters)
end

"""
    static_self_energy(swt::SpinWaveTheory, qpts, terms2)

Band-resolved energy shift caused by a correction `terms2` to the quadratic
Hamiltonian, as returned by [`hartree_fock_correction`](@ref) or
[`tadpole_correction`](@ref). Corrections from multiple sources should be
concatenated. The result has the same shape as [`dispersion`](@ref), to which it
is to be added.

Whereas [`corrected_dispersion`](@ref) rediagonalizes the corrected Hamiltonian,
this function evaluates the correction to first order only. The two differ at
relative order ``1/s²``, but only the latter can be combined with
[`cubic_self_energy`](@ref), which is likewise a first-order energy shift. It is
also the appropriate choice when the structure supports a Goldstone mode: an
``O(1/s)`` correction to a Hamiltonian with a protected zero mode produces a gap
of order ``\\sqrt{1/s}`` upon rediagonalization, obscuring the cancellation
between the terms that keeps the mode gapless.
"""
function static_self_energy(swt::SpinWaveTheory, qpts, terms2)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    H = zeros(ComplexF64, 2L, 2L)
    δH = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L)

    ret = zeros(L, length(qpts.qs))
    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        dynamical_matrix!(H, swt, q_reshaped)
        bogoliubov!(T, H)
        δH .= 0
        accum_quadratic!(δH, terms2, q_reshaped)
        # Writing the quadratic form as (1/2) y† (T† δH T) y in the quasi-particle
        # basis, the coefficient of α†_n α_n is the n-th diagonal element.
        view(ret, :, iq) .= real.(diag(T' * δH * T))[1:L]
    end

    return reshape(ret, L, size(qpts.qs)...)
end

"""
    corrected_dispersion(swt::SpinWaveTheory, qpts, terms2)

Excitation energies including a mean-field correction to the quadratic
Hamiltonian, as returned by [`hartree_fock_correction`](@ref) or
[`tadpole_correction`](@ref). Corrections from multiple sources should be
concatenated. Otherwise like [`dispersion`](@ref). See
[`static_self_energy`](@ref) for the alternative of applying the correction to
first order only.
"""
function corrected_dispersion(swt::SpinWaveTheory, qpts, terms2)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L)

    disp = zeros(L, length(qpts.qs))
    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        dynamical_matrix!(H, swt, q_reshaped)
        accum_quadratic!(H, terms2, q_reshaped)
        view(disp, :, iq) .= view(bogoliubov!(T, H), 1:L)
    end

    return reshape(disp, L, size(qpts.qs)...)
end
