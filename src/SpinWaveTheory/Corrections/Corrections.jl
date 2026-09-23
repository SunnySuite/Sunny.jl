# Corrections to linear spin wave theory (LSWT) at sub-leading order in 1/s.
#
# Expanding the spin Hamiltonian in Holstein-Primakoff bosons b, b† yields
#
#     H = E_cl + H₁ + H₂ + H₃ + H₄ + …,
#
# where Hₙ collects the n-boson terms and scales as s^(2-n/2). At a classical
# energy minimum H₁ vanishes. LSWT retains H₂ only, which `bogoliubov!`
# diagonalizes into quasi-particle modes α, α† with energies ε. Relative to
# LSWT, both H₄ (at first order) and H₃ (at second order) contribute at O(1/s).
#
# Writing quasi-particle momenta as subscripts and letting N denote the number
# of magnetic cells, the cubic term splits into a decay and a source channel,
#
#     H₃ = (1/2!√N) Σ_{1+2=3} [Γ₁(3; 1, 2) α†₃ α₁ α₂ + h.c.]
#        + (1/3!√N) Σ_{1+2+3=0} [Γ₂(1, 2, 3) α†₁ α†₂ α†₃ + h.c.],
#
# whose self-energies SelfEnergy.jl writes down. Only the decay channel has an
# imaginary part at T = 0, giving magnons a finite lifetime and transferring
# spectral weight into the two-magnon continuum.
#
# Sunny's Nambu conventions are reused throughout. In particular the columns
# `L+1:2L` of a Bogoliubov matrix `T` obtained at wavevector `q` are the
# eigenvectors at `-q` (see `excitations!`), so a band index ranging over the
# full Nambu space `1:2L` reaches both channels at once, which is how
# SelfEnergy.jl handles them together.
#
# Corrections are resummed rather than added, which conserves spectral weight.
# Collecting the quasi-particles into the Nambu vector y_𝐪 = [α_𝐪; α†_{-𝐪}]
# and writing Ĩ = diagm([ones(L), -ones(L)]) for the para-unitary metric, each
# correction enters the retarded Green function as a self-energy,
#
#     G(𝐪, ω) = (ω - diag(ε_𝐪) - Σ̂(𝐪, ω))⁻¹ Ĩ,
#
# where ε_𝐪 are the 2L signed energies returned by `bogoliubov!`, so that L
# poles at ω = ε_{𝐪n} are accompanied by L at ω = -ε_{-𝐪n}. The static mean
# fields of HartreeFock.jl, Tadpole.jl and `anisotropy_correction` contribute Σ̂
# = Ĩ (T†δH T)ᵗ for a perturbation (1/2) x†δH x of the quadratic Hamiltonian;
# SelfEnergy.jl contributes the frequency-dependent cubic self-energy. The
# transpose is required to make ĨΣ̂ Hermitian, without which the resummation
# would not preserve total weight. Both forms are verified against exact Green
# functions of a dimer.
#
# Observables are corrected too, by Observables.jl. Writing ũ[m, μ] for the
# amplitude with which observable μ creates Nambu mode m, obtained from the
# vectors u of `set_swt_observable_vectors!` as ũ = Tᵗ conj(u), the structure
# factor of CorrectedIntensities.jl is
#
#     S^{μν}(𝐪, ω) = -(1/π) Σ_{n,n′ ≤ L} ũ[n, μ] Im[G(𝐪, ω)][n, n′] conj(ũ[n′, ν]),
#
# where Im of a matrix means its anti-Hermitian part (G - G†)/2i. With Σ̂ = 0
# this reproduces the delta functions of `intensities_bands`, broadened.
#
# Only the block of G with n, n′ ≤ L appears, and it is obtained by projecting
# the Dyson equation onto that block, as in Eq. (12) of arXiv:1306.1231,
#
#     G_pp(𝐪, ω) = (ω - diag(ε_𝐪)_pp - Σ̂_pp(𝐪, ω))⁻¹,
#
# rather than by inverting the full 2L matrix and discarding the rest of the
# solution. The other blocks carry mirror poles at ω = -ε_{-𝐪n}, which a
# retarded spectral function weights negatively; at s = 1/2 a correction
# comparable to ε can push one up through ω = 0, making the full denominator
# near-singular, and the anomalous blocks of Σ̂ then carry a spurious pole into
# the particle block. Projecting also makes the result a spectral function in
# its own right: given Im Σ̂_pp ⪯ 0, which SelfEnergy.jl arranges, the
# denominator is nonsingular at every real ω, so S(𝐪, ω) ≥ 0 and is bounded by
# 1/πΓ, and since the denominator grows as ωI its frequency integral is the
# identity, conserving transverse weight exactly.
#
# The frequency-dependent momentum integrals are all one integral, over pairs of
# magnon lines at 𝐩 and 𝐪-𝐩, of the form
#
#     ∫d𝐤 Σ_{n₁n₂} V(𝐤, n₁, n₂) g(ω, x(𝐤, n₁, n₂)),
#
# with V ⪰ 0 independent of frequency and all the frequency dependence in a
# kernel g of the scalar pair energy x — a Cauchy denominator 1/(ω - x) in
# SelfEnergy.jl, a Lorentzian of half-width η in the two-magnon channel of
# CorrectedIntensities.jl. So V is accumulated into
# bins of x on the uniform grid of `loop_wavevectors` and g applied afterwards.
# That makes Im Σ̂_pp ⪯ 0 a property of the quadrature: `bin_index` splits each
# contribution between neighbouring bins with nonnegative weights, so a sum of
# positive semidefinite V stays positive semidefinite bin by bin, on any grid
# and at any s. Linear splitting preserves the zeroth and first moments exactly,
# so the sum rule survives binning; the shape error is O((Δ/Γ)²).
#
# It is one integral rather than two because a two-magnon final state is
# reachable by two interfering routes: the transverse (odd) part of the
# observable creates one magnon which the cubic vertex splits, while the
# longitudinal part Sᶻ = s - b†b creates the pair directly; see Observables.jl.
# Both are of order s⁰ relative to the one-magnon amplitude, so the interference
# is of the same relative order 1/s as either squared. Writing v for the vertex
# factor and β for the pair amplitude, V is the rank-one yy† built from
#
#     y = [√18 v ; β],
#
# whose diagonal blocks are the cubic self-energy and the two-magnon continuum
# and whose off-diagonal block is the interference. This interference term was
# largely omitted from the relevant literature, e.g., the sequence of works by
# Chernyshev et al. (arXiv:1306.1231 and arXiv:1607.08238). In that context the
# interference effects are tiny -- a couple-percent error in the two-magnon
# weight, which is itself a small fraction of the total intensity. The reason is
# a sum rule, Σ_𝐪 ∫dω cross = 0, which holds for a trace measure because 𝐒·𝐒
# links no odd number of bosons to an even one. It makes the interference cancel
# rather than making it small: pointwise it is comparable to the continuum, and
# only the integral collapses. A readout that is not a trace has no such
# protection, and there the interference can be huge, e.g., for a chiral readout
# such as Im S^{xy} of the XXZ model in an out-of-plane field. There it exceeds
# the two-magnon continuum it interferes with, and the direct channel is five
# orders of magnitude below either.
#
# The test suite certifies all the correction terms above by comparing to exact
# calculations on a dimer model with arbitrary anisotropic interactions and
# readouts.

# Why the 1/s corrections of this directory are unavailable for `swt`, or `nothing`
# if they are available. Returned rather than thrown so that a caller offering a
# correct but weaker result in the unsupported cases, such as
# `corrected_magnetic_moments`, can ask without catching.
function corrections_unsupported_reason(swt::SpinWaveTheory)
    (; sys) = swt
    @assert sys.mode in (:dipole, :dipole_uncorrected, :SUN)

    sys.mode == :SUN && return "are not yet implemented in :SUN mode"
    is_entangled(sys) && return "are not supported for entangled units"
    isnothing(sys.ewald) || return "do not yet support long-range dipole-dipole interactions"

    for int in sys.interactions_union
        for pc in int.pair
            pc.isculled && break
            iszero(pc.biquad) || return "do not yet support biquadratic exchange"
        end
    end
    return nothing
end

# Errors unless `swt` describes a model for which the 1/s corrections in this
# directory are implemented.
function check_corrections_supported(swt::SpinWaveTheory)
    reason = corrections_unsupported_reason(swt)
    isnothing(reason) || error("1/s corrections $reason.")
end

# Wavevectors 𝐩 of the loop integrals over the magnetic Brillouin zone, for an
# integrand that pairs a line at 𝐩 with one at 𝐪-𝐩. Each channel diverges at
# the zone centre, so the grid must avoid it in 𝐩 and in 𝐪-𝐩 alike;
# offsetting by half a step does only the former. In units of a step, and per
# dimension, the forbidden offsets are 0, putting 𝐩 on the zone centre, and t =
# dims*𝐪 mod 1, putting 𝐪-𝐩 there. Sit at the midpoint of the larger arc
# between them, which keeps a quarter step of clearance on both lines and
# reduces to the half step when 𝐪 is commensurate with the grid. Without this
# the integral picks up a spurious divergence from one grid point whenever
# dims*𝐪 has a half-integer component, afflicting isolated wavevectors of a
# path rather than all of them.
function loop_wavevectors(dims, q_reshaped=zero(Vec3))
    offsets = ntuple(3) do d
        t = mod(dims[d] * q_reshaped[d], 1)
        t < 1/2 ? (t + 1)/2 : t/2
    end
    return [Vec3((i - 1 + offsets[1]) / dims[1], (j - 1 + offsets[2]) / dims[2], (k - 1 + offsets[3]) / dims[3])
            for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]]
end

# The wavevector loop shared by every frequency-dependent momentum integral of
# this module, as described above. Calls `f(𝐩, T1, T2, ε1, ε2)` once per
# wavevector 𝐩 of `ps`, where T1 = T(𝐩) and T2 = T(𝐪-𝐩) diagonalize the two
# magnon lines being paired, and ε1, ε2 are the signed energies that
# `bogoliubov!` returns with them. All of these are overwritten on each
# iteration, so `f` must consume them before returning.
function foreach_magnon_pair(f, swt::SpinWaveTheory, q_reshaped, ps)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    T1 = zeros(ComplexF64, 2L, 2L)
    T2 = zeros(ComplexF64, 2L, 2L)
    for p in ps
        dynamical_matrix!(H, swt, p)
        ε1 = bogoliubov!(T1, H)
        dynamical_matrix!(H, swt, q_reshaped - p)
        ε2 = bogoliubov!(T2, H)
        f(p, T1, T2, ε1, ε2)
    end
end

# Amplitude for the longitudinal part of observable μ to create the pair of
# magnons (𝐩 a, 𝐪-𝐩 b) directly, via Sᶻ = s - b†b. Only the component of each
# observable along the local quantization axis contributes, the transverse ones
# being odd in the boson number; `pref[μ, i]` carries that component together
# with the phase factor of site i, as `observable_prefactor` defines it. The two
# terms symmetrize over which line takes the b† of b†b, and the 1/√2 is the norm
# of the symmetrized two-boson state.
#
# Both lines are expressed in the Bogoliubov matrices at +𝐩 and +(𝐪-𝐩), which
# `foreach_magnon_pair` supplies and the cubic vertex shares. An equivalent form
# in T(-𝐩) and T(𝐩-𝐪) follows from the Nambu symmetry of SelfEnergy.jl, but
# would come from independent diagonalizations, whose free per-band phase the
# interference cannot tolerate. The leading minus sign, that of Sᶻ = s - b†b,
# cancels in the |β|² of the direct channel but is the whole sign of the
# interference.
function pair_amplitude(pref, T1, T2, a, b, μ, L)
    return -sum(1:L) do i
        pref[μ, i] * (T1[i, a]*T2[L+i, b] + T1[L+i, a]*T2[i, b])
    end / √2
end

# Longitudinal observable prefactors for `pair_amplitude`, at one wavevector.
function pair_amplitude_prefactors!(pref, swt::SpinWaveTheory, q_reshaped, q_global)
    (; sys, measure, data) = swt
    for μ in 1:num_observables(measure), i in 1:nbands(swt)
        O = (data::SWTDataDipole).observables[μ, i]
        pref[μ, i] = conj(observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)) * O[3]
    end
    return pref
end

# Applies `f` to each index, optionally in parallel. The wavevector loops of this
# module allocate their buffers per iteration so that they may be threaded.
function foreach_maybe_threaded(f, threaded, indices)
    if threaded
        Threads.@threads for i in indices
            f(i)
        end
    else
        foreach(f, indices)
    end
end

# Dimensions of the loop grid needed to reach a relative accuracy `tol` at
# regulator `η`. Every frequency-dependent integrand here is a function of the
# pair energy x(𝐤) = ε_𝐤 + ε_{𝐪-𝐤} smoothed on the scale η, so the grid must
# resolve x to within η. The number of points along a direction therefore goes
# as the range x sweeps there divided by η, estimated below by the range each
# band sweeps along a line, doubled for the two magnons; a non-dispersing
# direction needs no grid.
#
# The dependence on `tol` and the prefactor are calibrated rather than derived,
# convergence being algebraic because the dispersion is non-analytic at the
# Goldstone wavevectors. Measured on the triangular-lattice antiferromagnet, the
# error in the integrated weight falls as n^-1.6 with a factor-of-two scatter,
# since how closely the grid approaches a near-singular point depends on n
# arithmetically. The prefactor carries margin accordingly. The cost grows as
# 1/√tol per dimension, so a tenfold tighter tolerance is a tenfold longer
# calculation in two dimensions.
function auto_loop_grid(swt::SpinWaveTheory, η, tol)
    ncoarse = 8
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L)
    ε = zeros(L, ncoarse, ncoarse, ncoarse)
    for i in 1:ncoarse, j in 1:ncoarse, k in 1:ncoarse
        q = Vec3((i - 1/2)/ncoarse, (j - 1/2)/ncoarse, (k - 1/2)/ncoarse)
        dynamical_matrix!(H, swt, q)
        view(ε, :, i, j, k) .= view(bogoliubov!(T, H), 1:L)
    end

    return ntuple(3) do d
        r = maximum(maximum(ε; dims=d+1) - minimum(ε; dims=d+1))
        # The denominator is 0.8η at the default tolerance, and shrinks as √tol
        2r < η ? 1 : max(4, ceil(Int, 2r / (8η * √tol)))
    end
end

# Index `b` and interpolation weight `f` for scattering a pair energy `x ≥ 0`
# into bins centered at (b - 1)Δ, b = 1, 2, …: a mass m at x becomes (1-f)m in
# bin b and f m in bin b+1. See the discussion of binning above. The accumulator
# `ρ` is grown to length b+1 with new bins made by `mk`, so that the caller need
# only add into them.
function bin_index!(ρ, x, Δ, mk)
    t = max(x, 0) / Δ
    b = 1 + floor(Int, t)
    while length(ρ) < b + 1
        push!(ρ, mk())
    end
    return (b, t - (b - 1))
end
