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
# giving self-energies
#
#     Σᵈ(k, ω) = (1/2N) Σ_q Σ_{n₁n₂} |Γ₁(k; q n₁, k-q n₂)|² / (ω - ε_{q n₁} - ε_{k-q n₂} + iη)
#     Σˢ(k, ω) = -(1/2N) Σ_q Σ_{n₁n₂} |Γ₂(k, -q n₁, q-k n₂)|² / (ω + ε_{q n₁} + ε_{k-q n₂}).
#
# Only Σᵈ has an imaginary part at T = 0; it is what gives magnons a finite
# lifetime and transfers spectral weight into the two-magnon continuum.
#
# Sunny's Nambu conventions are reused throughout. In particular the columns
# `L+1:2L` of a Bogoliubov matrix `T` obtained at wavevector `q` are the
# eigenvectors at `-q` (see `excitations!`), so a band index ranging over the
# full Nambu space `1:2L` reaches both channels above.
#
# Corrections are resummed rather than added, which is what conserves spectral
# weight. Collecting the quasi-particles into the Nambu vector y_𝐪 = [α_𝐪; α†_{-𝐪}]
# and writing Ĩ = diagm([ones(L), -ones(L)]) for the para-unitary metric, the
# retarded Green function of LSWT is
#
#     G₀(𝐪, ω) = (ω - diag(ε_𝐪))⁻¹ Ĩ,
#
# where ε_𝐪 are the 2L signed energies returned by `bogoliubov!`, so that the L
# poles at ω = ε_{𝐪n} are accompanied by L poles at ω = -ε_{-𝐪n}. Each correction
# enters as a self-energy,
#
#     G(𝐪, ω) = (ω - diag(ε_𝐪) - Σ̂(𝐪, ω))⁻¹ Ĩ,
#
# with contributions from the static mean fields of HartreeFock.jl, Tadpole.jl and
# `anisotropy_correction`, which for a perturbation (1/2) x†δH x of the quadratic
# Hamiltonian take the form Σ̂ = Ĩ (T†δH T)ᵗ, and from the frequency-dependent
# cubic self-energy of SelfEnergy.jl. The transpose is not cosmetic; it is what
# makes ĨΣ̂ Hermitian, which in turn is what makes the resummation preserve total
# weight. Both forms are verified against exact Green functions of a dimer.
#
# Observables are corrected too, by Observables.jl. Writing ũ[m, μ] for the
# amplitude with which observable μ creates Nambu mode m, obtained from the
# vectors u of `set_swt_observable_vectors!` as ũ = Tᵗ conj(u), the structure
# factor of CorrectedIntensities.jl is
#
#     S^{μν}(𝐪, ω) = -(1/π) Σ_{n,n′ ≤ L} ũ[n, μ] Im[G(𝐪, ω)][n, n′] conj(ũ[n′, ν]),
#
# where Im of a matrix means its anti-Hermitian part (G - G†)/2i. With Σ̂ = 0 this
# reproduces the delta functions of `intensities_bands`, broadened.
#
# Only the block of G with n, n′ ≤ L appears. The other blocks carry the poles at
# ω = -ε_{-𝐪n}, which the spectral function of a retarded correlator weights
# negatively; they are the mirror images of the physical excitations, present so
# that S(𝐪, ω) may be continued to ω < 0, and including them would subtract
# weight at ω > 0 through the tails of the resolution function.
#
# That block is obtained by projecting the Dyson equation onto n, n′ ≤ L rather than
# by inverting the full 2L matrix and discarding the rest of the solution,
#
#     G_pp(𝐪, ω) = (ω - diag(ε_𝐪)_pp - Σ̂_pp(𝐪, ω))⁻¹,
#
# as in Eq. (12) of Mourigal et al., PRB 88, 094407 (2013). The two agree to the
# order worked to, since a mirror pole can reach the particle block only through the
# anomalous blocks of Σ̂ twice over, at O(1/s²). They do not agree at s = 1/2, where a
# correction comparable to ε can push a mirror pole up through ω = 0 to collide with a
# physical one; the pair then leaves the real axis, the full 2L denominator becomes
# near-singular, and the anomalous blocks of Σ̂ carry the resulting spurious pole into
# the particle block.
#
# Projecting the equation is also what makes the corrected structure factor a
# spectral function in its own right, rather than one only to the order worked to.
# Provided Im Σ̂_pp ⪯ 0, which SelfEnergy.jl arranges, the denominator has imaginary
# part ⪰ ΓI and so is nonsingular at every real ω, making S(𝐪, ω) non-negative and
# bounded by the resolution height 1/πΓ. And since the denominator grows as ωI, the
# frequency integral of the spectral function of G_pp is the identity, so the magnon
# poles may move and broaden but Σ_n |ũ[n, μ]|² of Observables.jl remains exactly the
# whole transverse weight, which is what the sum rule requires.
#
# Two of the momentum integrals here depend on frequency: the cubic self-energy of
# SelfEnergy.jl and the two-magnon continuum of TwoMagnon.jl. Both have the form
#
#     ∫d𝐤 Σ_{n₁n₂} V(𝐤, n₁, n₂) g(ω, x(𝐤, n₁, n₂)),
#
# with V ⪰ 0 independent of frequency and all the frequency dependence in a kernel
# g of the single scalar pair energy x — a Cauchy denominator 1/(ω - x) in the first
# case, the resolution kernel in the second. Both are therefore evaluated by
# accumulating V into bins of x on the uniform grid of `loop_wavevectors`, and applying
# g afterwards. The wavevector loop then costs nothing per frequency, and it is what
# makes Im Σ̂_pp ⪯ 0 a property of the quadrature rather than an accident of it:
# `bin_index` splits each contribution between two neighbouring bins with weights
# that are nonnegative by construction, so a sum of positive semidefinite V stays
# positive semidefinite bin by bin, on any grid and at any s.
#
# The price is a discretization of x. Splitting a mass linearly between neighbours
# preserves both the zeroth and the first moment of the binned measure exactly, so
# the large-ω tail Σ̂ → (Σ_b ρ_b)/ω is untouched and the ∫dω A = I sum rule above
# survives binning identically; the leading error is O(Δ²) times the curvature of g,
# hence O((Δ/Γ)²), which a bin width Δ = fwhm/32 puts some fifty times below the
# error of the wavevector grid itself.

# Errors unless `swt` describes a model for which the 1/s corrections in this
# directory are implemented.
function check_corrections_supported(swt::SpinWaveTheory)
    (; sys) = swt

    sys.mode == :SUN && error("1/s corrections are not yet implemented in :SUN mode.")
    @assert sys.mode in (:dipole, :dipole_uncorrected)
    is_entangled(sys) && error("1/s corrections are not supported for entangled units.")
    isnothing(sys.ewald) || error("1/s corrections do not yet support long-range dipole-dipole interactions.")

    for int in sys.interactions_union
        for pc in int.pair
            pc.isculled && break
            iszero(pc.biquad) || error("1/s corrections do not yet support biquadratic exchange.")
        end
    end
end

# Wavevectors 𝐩 of the loop integrals over the magnetic Brillouin zone, for an
# integrand that pairs a line at 𝐩 with one at 𝐪-𝐩. Both lines are singular at the
# zone centre, the Goldstone wavevector of the magnetic cell, where the integrand
# itself is finite but each of its channels diverges, so the grid must avoid that
# point in 𝐩 and in 𝐪-𝐩 alike. Offsetting by half a step does only the former. In
# units of a step, and per dimension, the forbidden offsets are 0, which puts 𝐩 on the
# zone centre, and t = dims*𝐪 mod 1, which puts 𝐪-𝐩 there. Sit at the midpoint of the
# larger of the two arcs between them, which keeps a quarter step of clearance on both
# lines and reduces to the half step when 𝐪 is commensurate with the grid.
#
# Without this the loop integral develops a spurious divergent contribution from a
# single grid point whenever dims*𝐪 has a half-integer component. It is easy to miss,
# because it afflicts isolated wavevectors of a path rather than all of them.
function loop_wavevectors(dims, q_reshaped=zero(Vec3))
    offsets = ntuple(3) do d
        t = mod(dims[d] * q_reshaped[d], 1)
        t < 1/2 ? (t + 1)/2 : t/2
    end
    return [Vec3((i - 1 + offsets[1]) / dims[1], (j - 1 + offsets[2]) / dims[2], (k - 1 + offsets[3]) / dims[3])
            for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]]
end

# Dimensions of the loop grid needed to reach a relative accuracy `tol` at regulator
# `η`. Every frequency-dependent integrand here is a function of the pair energy
# x(𝐤) = ε_𝐤 + ε_{𝐪-𝐤} smoothed on the scale η, so what the grid must do is resolve x
# to within η. That fixes the scaling: the number of points along a direction goes as
# the range that x sweeps along that direction divided by η, estimated below by the
# range that each band sweeps along a line, doubled because the pair energy involves
# two magnons. A direction along which the magnons do not disperse needs no grid.
#
# The dependence on `tol`, and the prefactor, are calibrated rather than derived.
# Convergence is algebraic and not exponential, because the dispersion is itself
# non-analytic at the Goldstone wavevectors of an ordered structure. Measured on the
# triangular-lattice antiferromagnet, the error in the integrated weight falls as
# n^-1.6, but with a factor-of-two scatter about that trend, since how closely the grid
# approaches a near-singular point depends on n arithmetically rather than smoothly.
# The prefactor therefore carries margin: at the default tol = 0.01 the grid comes out a
# little finer than the one hand-tuned for Fig. 2 of Mourigal et al., and the achieved
# error in the integrated weight is around half of `tol`. Note what tightening costs.
# The number of points grows as 1/√tol per dimension, so a tenfold tighter tolerance is
# a tenfold longer calculation in two dimensions.
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

# Index `b` and interpolation weight `f` for scattering a pair energy `x ≥ 0` into
# bins centered at (b - 1)Δ, b = 1, 2, …: a mass m at x becomes (1-f)m in bin b and
# f m in bin b+1. See the discussion of binning above; the caller grows its own
# accumulator to length b+1.
function bin_index(x, Δ)
    t = max(x, 0) / Δ
    b = 1 + floor(Int, t)
    return (b, t - (b - 1))
end
