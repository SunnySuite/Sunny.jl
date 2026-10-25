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
# whose self-energies SelfEnergy.jl writes down. At T = 0 only the decay channel
# resonates at ω > 0, giving magnons a finite lifetime.
#
# Sunny's Nambu conventions are reused throughout. In particular the columns
# `L+1:2L` of a Bogoliubov matrix `T` obtained at wavevector `q` are the
# eigenvectors at `-q` (see `excitations!`), so a band index ranging over the
# full Nambu space `1:2L` reaches both channels at once.
#
# Resummation. Collect the quasi-particles into y_𝐪 = [α_𝐪; α†_{-𝐪}], with
# metric Ĩ = diagm([ones(L); -ones(L)]), and write w = T†u for the observable
# amplitudes (u from `set_swt_observable_vectors!`, corrected by
# Observables.jl). LSWT gives the retarded response
#
#     χ^{μν}(z) = w_μ' G₀(z) w_ν,   G₀(z) = (zĨ - |ε|)⁻¹,
#
# whose anti-Hermitian part (χ' - χ)/2πi is the broadened `intensities` at z = ω
# + iη, including the mirror poles at ω < 0.
#
# At O(1/s) the static mean fields of HartreeFock.jl, Tadpole.jl and
# `anisotropy_correction` add Σstat = T†δH T, for a perturbation (1/2)x†δH x of
# the quadratic Hamiltonian. The cubic vertex couples each magnon to pairs of
# magnons, and the observable creates pairs directly too (Sᶻ = s - b†b). Both
# are captured by an auxiliary quadratic model: magnons coupled to a bath of
# free two-magnon states. For each pair of internal lines (𝐩 a, 𝐪-𝐩 b) at
# pair energy x, define
#
#     y = [√18 U[a, b, :]; β],
#
# the vertex to each of the 2L external Nambu legs and the amplitudes β for each
# observable to create the pair (`pair_amplitude`). Forward lines (x > 0) are
# bath particles, backward lines (x < 0) bath holes of the opposite metric sign.
# Integrating out the bath gives the Cauchy transform
#
#     K(z) = Σ_pairs ± y y† / (z - x),
#
# with blocks K_mm (the cubic self-energy), K_md, K_dm and K_dd (the direct
# two-magnon continuum). The exact response of the auxiliary model is then
#
#     χ = w'Gw + K_dm G w + w'G K_md + K_dm G K_md + K_dd,
#     G = (zĨ - |ε| - Σstat - K_mm)⁻¹.
#
# This is exact at O(1/s) and, being the resolvent of a quadratic model,
# inherits its structure: Nambu symmetry, Goldstone protection (the static and
# dynamic 1/ε divergences cancel in the full 2L inverse), η as a pure
# Lorentzian, and the commutator sum rule ∫dω S = w'Ĩw. It is positive at ω > 0
# whenever the auxiliary model is stable, i.e. |ε| + Σstat + K_mm(0) ⪰ 0. Where
# it is not, the 1/s correction to some mode is comparable to its harmonic
# energy, and a resummed pole moves onto the imaginary axis. This is a breakdown
# of the expansion for that mode rather than a physical instability, and that
# mode may carry little weight in the observable; `corrected_channels` reports
# it. The interference terms linear in K_md are absent from the standard 1/s
# treatment, which takes the magnon spectral function as the major component of
# S and adds the two-magnon continuum separately (PRB 79, 144416, Sec. VI). For
# a trace measure they cancel in a zone sum, but not pointwise.
#
# K depends on frequency only through the scalar x, so the masses y y†, summed
# over the wavevectors of `loop_wavevectors`, are accumulated into bins of x and
# the Cauchy transform is applied afterwards. Linear splitting between
# neighbouring bins keeps each channel's measure semidefinite and preserves its
# zeroth and first moments; the shape error is O((Δ/η)²).
#
# The test suite certifies every term by comparing to exact diagonalization of a
# cluster with anisotropic interactions and readouts.

# Why the 1/s corrections of this directory are unavailable for `swt`, or `nothing`
# if they are available. Returned rather than thrown so that a caller offering a
# correct but weaker result in the unsupported cases, such as
# `corrected_magnetic_moments`, can ask without catching.
function corrections_unsupported_reason(swt::SpinWaveTheory)
    (; sys) = swt
    @assert sys.mode in (:dipole, :dipole_uncorrected, :SUN)

    is_entangled(sys) && return "are not supported for entangled units"
    isnothing(sys.ewald) || return "do not yet support long-range dipole-dipole interactions"

    # In :SUN mode every coupling, biquadratic included, has been decomposed into
    # the tensor pairs that `sun_monomials` expands, so there is nothing to reject.
    if sys.mode != :SUN
        for int in sys.interactions_union
            for pc in int.pair
                pc.isculled && break
                iszero(pc.biquad) || return "do not yet support biquadratic exchange"
            end
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
# integrand that pairs a line at 𝐩 with one at 𝐪-𝐩, carrying a multiplicity
# `wts` each and totalling `npts = prod(dims)` points. Each channel diverges at
# the zone centre, so the grid must avoid it in 𝐩 and in 𝐪-𝐩 alike;
# offsetting by half a step does only the former. In units of a step, and per
# dimension, the forbidden offsets are 0, putting 𝐩 on the zone centre, and t =
# dims*𝐪 mod 1, putting 𝐪-𝐩 there. Sit at the midpoint of the larger arc
# between them, which keeps a quarter step of clearance on both lines and
# reduces to the half step when 𝐪 is commensurate with the grid. Without this
# the integral picks up a spurious divergence from one grid point whenever
# dims*𝐪 has a half-integer component, afflicting isolated wavevectors of a
# path rather than all of them.
#
# That offset also makes the grid closed under the involution 𝐩 ↦ 𝐪-𝐩, which
# halves the work: the two magnon lines are interchangeable, so a point and its
# partner contribute equally to every integrand of this module, and only one of
# the two need be visited. Closure is exact rather than approximate, since
# dims*𝐪 - 2*offset is dims*𝐪 minus its own fractional part, less one in the
# branch that adds a half step. Writing 𝐦 for that integer, the partner of grid
# index 𝐢 is 𝐦 - 𝐢 + 2 taken mod dims; the fixed points of the map, at most one
# per dimension pair, keep multiplicity one.
struct LoopGrid
    ps::Vector{Vec3}
    wts::Vector{Float64}
    npts::Int
end

function loop_wavevectors(dims, q_reshaped=zero(Vec3))
    offsets = ntuple(3) do d
        t = mod(dims[d] * q_reshaped[d], 1)
        t < 1/2 ? (t + 1)/2 : t/2
    end
    ms = ntuple(d -> round(Int, dims[d] * q_reshaped[d] - 2offsets[d]), 3)
    partner = c -> CartesianIndex(ntuple(d -> mod(ms[d] - c[d] + 1, dims[d]) + 1, 3))

    (ps, wts) = (Vec3[], Float64[])
    visited = falses(dims)
    for c in CartesianIndices(visited)
        visited[c] && continue
        visited[c] = visited[partner(c)] = true
        push!(ps, Vec3(ntuple(d -> (c[d] - 1 + offsets[d]) / dims[d], 3)))
        push!(wts, partner(c) == c ? 1 : 2)
    end
    return LoopGrid(ps, wts, prod(dims))
end

# The wavevector loop shared by every frequency-dependent momentum integral of
# this module, as described above. Calls `f(𝐩, w, T1, T2, ε1, ε2)` once per
# wavevector 𝐩 of `grid`, where T1 = T(𝐩) and T2 = T(𝐪-𝐩) diagonalize the two
# magnon lines being paired, ε1, ε2 are the signed energies that `bogoliubov!`
# returns with them, and `w` is the multiplicity by which the contribution is to
# be scaled. All of the arrays are overwritten on each iteration, so `f` must
# consume them before returning.
function foreach_magnon_pair(f, swt::SpinWaveTheory, q_reshaped, grid::LoopGrid)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    # Both magnon lines are live at once, so each gets its own workspace
    (ws1, ws2) = (BogoliubovWorkspace(L), BogoliubovWorkspace(L))
    for (p, w) in zip(grid.ps, grid.wts)
        dynamical_matrix!(H, swt, p)
        ε1 = bogoliubov!(ws1, H)
        dynamical_matrix!(H, swt, q_reshaped - p)
        ε2 = bogoliubov!(ws2, H)
        f(p, w, ws1.T, ws2.T, ε1, ε2)
    end
end

# Site owning boson `a` of the Nambu labeling, bosons being laid out as (flavor,
# atom) with flavor fastest, `Nf` of them per site. The identity in dipole mode.
boson_site(a, L, Nf) = div(mod1(a, L) - 1, Nf) + 1

# Amplitude for the even part of an observable to create the pair of magnons (𝐩
# a, 𝐪-𝐩 b) directly, given that part as the `BosonMonomial{2}` list that
# `observable_pair_words` builds, already carrying its Fourier phase. That type
# is declared in Vertices.jl, which is included after this file, so it cannot
# appear in the signature. Only the even part of an observable contributes here,
# the odd one changing the boson number by one. The two terms symmetrize over
# which line takes which slot of the word, and the 1/√2 is the norm of the
# symmetrized two-boson state. No phase accompanies the slots because an
# observable is onsite, so both of them carry the same cell offset.
#
# Both lines are expressed in the Bogoliubov matrices at +𝐩 and +(𝐪-𝐩), which
# `foreach_magnon_pair` supplies and the cubic vertex shares. An equivalent form
# in T(-𝐩) and T(𝐩-𝐪) follows from the Nambu symmetry of SelfEnergy.jl, but
# would come from independent diagonalizations, whose free per-band phase the
# interference cannot tolerate. The sign of the word, that of Sᶻ = s - b†b in
# dipole mode, cancels in the |β|² of the direct channel but is the whole sign
# of the interference.
function pair_amplitude(words, T1, T2, a, b)
    return sum(words; init=zero(ComplexF64)) do (; c, as)
        c * (T1[as[2], a]*T2[as[1], b] + T1[as[1], a]*T2[as[2], b])
    end / √2
end

# Applies `f` to each index, optionally in parallel. The wavevector loops of
# this module allocate their buffers per iteration so that they may be threaded.
# A progress bar labeled `desc` is shown unless `desc` is nothing; `next!` is
# itself thread safe. To animate the bar, stdout must allow the `\r` character
# to rewrite the current line; this is only supported on TTY outputs.
function foreach_maybe_threaded(f, threaded, indices; desc=nothing)
    enabled = !isnothing(desc) && stdout isa Base.TTY
    meter = ProgressMeter.Progress(length(indices); desc=@something(desc, ""),
                                  enabled, output=stdout)
    g = i -> (f(i); ProgressMeter.next!(meter))
    if threaded
        Threads.@threads for i in indices
            g(i)
        end
    else
        foreach(g, indices)
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

# Matrix-valued measure over a bath energy x of either sign, binned as described
# above: bin j is centered at jΔ. Each channel gets its own measure. Bins are
# Hermitian, so each stores only its upper triangle, packed column by column
# (the BLAS "packed" layout), and is read with `bin_entry`.
struct PairMeasure
    Δ::Float64
    dim::Int
    bins::Dict{Int, Vector{ComplexF64}}
end

PairMeasure(Δ, dim) = PairMeasure(Δ, dim, Dict{Int, Vector{ComplexF64}}())

# Position of element (n, n′), n ≤ n′, in a packed upper triangle
packed_index(n, n′) = n + n′ * (n′ - 1) ÷ 2

# Element (n, n′) of packed Hermitian bin `v`, for any n, n′
bin_entry(v, n, n′) = n ≤ n′ ? v[packed_index(n, n′)] : conj(v[packed_index(n′, n)])

# Accumulates the rank-one mass c y y† at bath energy x. The weight c carries the
# sign of the channel, and would carry its thermal factor at T > 0.
function accum_binned!(ρ::PairMeasure, x, c, y)
    t = x / ρ.Δ
    j = floor(Int, t)
    f = t - j
    for (jj, cc) in ((j, c * (1 - f)), (j + 1, c * f))
        v = get!(() -> zeros(ComplexF64, packed_index(ρ.dim, ρ.dim)), ρ.bins, jj)
        @inbounds for n′ in 1:ρ.dim, n in 1:n′
            v[packed_index(n, n′)] += cc * y[n] * conj(y[n′])
        end
    end
end

# Cauchy transform K(z) = Σ_j ρ_j / (z - jΔ), summed over the measures `ρs` and
# returned as a `dim×dim×length(zs)` array. Frequencies are processed in small
# blocks, which keeps the matrix products cache resident.
function cauchy_transform(ρs, zs; nb=16)
    dim = first(ρs).dim
    K = zeros(ComplexF64, dim, dim, length(zs))
    Kr = reshape(K, dim^2, length(zs))
    for ρ in ρs
        js = collect(keys(ρ.bins))
        P = zeros(ComplexF64, dim^2, length(js))
        for (i, j) in enumerate(js), n′ in 1:dim, n in 1:dim
            P[n + (n′ - 1) * dim, i] = bin_entry(ρ.bins[j], n, n′)
        end
        C = zeros(ComplexF64, length(js), nb)
        for r in Iterators.partition(eachindex(zs), nb)
            Cr = view(C, :, 1:length(r))
            for (k, iz) in enumerate(r), (i, j) in enumerate(js)
                Cr[i, k] = 1 / (zs[iz] - j * ρ.Δ)
            end
            mul!(view(Kr, :, r), P, Cr, true, true)
        end
    end
    return K
end
