# Momentum-space quadrature: the Brillouin-zone averages of the static
# corrections, the loop grid of the frequency-dependent ones, and the binned
# measure over pair energy into which a loop is accumulated.

# ---- Brillouin-zone quadrature ----

"""
    BZGrid(n1, n2, n3)

A uniform grid of `n1×n2×n3` wavevectors in the magnetic Brillouin zone, offset
by half a step from the zone centre, for use in place of a numeric `tol` by the
1/s corrections, e.g. [`corrected_intensities`](@ref).

Every 1/s correction is a momentum integral, and `tol` selects how it is done.
A number is a relative accuracy target: static averages, e.g. mean fields and
energies, use adaptive cubature, and the frequency-dependent loop integrals use
a uniform grid that is fine enough to resolve the regulator `η`. A `BZGrid`
instead fixes every integral, static and loop alike, to an average over these
wavevectors. A shared grid makes identities between static and loop terms exact
grid by grid, e.g. the Ward identity that keeps a Goldstone mode gapless, and
it gives results that are smooth in the model parameters, as fitting requires.
"""
struct BZGrid
    dims :: NTuple{3, Int}
end

BZGrid(n1, n2, n3) = BZGrid((n1, n2, n3))

# Scale of the error that the quadrature leaves in an identity that holds only
# for exact integrals, against which such identities are asserted.
quadrature_noise(tol::Real) = max(tol, 1e-8)
quadrature_noise(::BZGrid) = 1e-3

# Wavevectors of a uniform grid, offset by half a step to keep off the zone
# centre, in reshaped RLU. The grid is closed under 𝐪 ↦ -𝐪, which keeps every
# average of a Hermitian quantity Hermitian.
grid_points(grid::BZGrid) = [Vec3((Tuple(c) .- 1/2) ./ grid.dims) for c in CartesianIndices(grid.dims)]

# Average of `f(q_reshaped)`, which may be array valued, over the magnetic
# Brillouin zone 𝐪 ∈ [0, 1)³ in reshaped RLU.
function bz_average(f, grid::BZGrid)
    ps = grid_points(grid)
    return sum(f, ps) / length(ps)
end

# HCubature stops once `err ≤ max(atol, rtol * norm(val))`. Setting `atol`
# equal to `tol` measures the accuracy against max(norm(val), 1), so that
# averages vanishing by symmetry converge at once. The averages taken here,
# correlations and energies per site, are of order one.
bz_average(f, tol::Real) = first(hcubature(q -> f(Vec3(q)), (0, 0, 0), (1, 1, 1); rtol=tol, atol=tol))

# The `tol` with which `vac` is to be used: its own, if it was solved on one,
# otherwise `tol`, falling back to `default`
function resolve_tol(swt, vac, tol, default=nothing)
    vac.swt === swt || error("Vacuum must be built on the same `SpinWaveTheory`")
    if !isnothing(vac.tol)
        isnothing(tol) || tol == vac.tol ||
            error("Vacuum was solved with `tol = $(vac.tol)`, which `tol = $tol` would contradict")
        return vac.tol
    end
    tol = @something tol default error("Must specify `tol` to control momentum-space integration; see `BZGrid`.")
    tol isa Union{Real, BZGrid} || error("`tol` must be a number or a `BZGrid`")
    return tol
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

# Dimensions of the loop grid needed to reach a relative accuracy `tol` at
# regulator `η`, for the internal lines of the vacuum `vac`. Every
# frequency-dependent integrand here is a function of the pair energy x(𝐤) =
# ε_𝐤 + ε_{𝐪-𝐤} smoothed on the scale η, so the grid must resolve x to within
# η. The number of points along a direction therefore goes as the range x sweeps
# there divided by η, estimated below by the range each band sweeps along a
# line, doubled for the two magnons; a non-dispersing direction needs no grid.
#
# The dependence on `tol` and the prefactor are calibrated rather than derived,
# convergence being algebraic because the dispersion is non-analytic at the
# Goldstone wavevectors. Measured on the triangular-lattice antiferromagnet, the
# error in the integrated weight falls as n^-1.6 with a factor-of-two scatter,
# since how closely the grid approaches a near-singular point depends on n
# arithmetically. The prefactor carries margin accordingly. The cost grows as
# 1/√tol per dimension, so a tenfold tighter tolerance is a tenfold longer
# calculation in two dimensions.
function auto_loop_grid(vac, η, tol)
    ncoarse = 8
    L = nbands(vac.swt)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    ε = zeros(L, ncoarse, ncoarse, ncoarse)
    for i in 1:ncoarse, j in 1:ncoarse, k in 1:ncoarse
        q = Vec3((i - 1/2)/ncoarse, (j - 1/2)/ncoarse, (k - 1/2)/ncoarse)
        view(ε, :, i, j, k) .= view(vacuum_bogoliubov!(ws, H, vac, q), 1:L)
    end

    return ntuple(3) do d
        r = maximum(maximum(ε; dims=d+1) - minimum(ε; dims=d+1))
        # The denominator is 0.8η at the default tolerance, and shrinks as √tol
        2r < η ? 1 : max(4, ceil(Int, 2r / (8η * √tol)))
    end
end

# Dimensions of the loop grid that `tol` selects: those of a `BZGrid`, or for a
# relative accuracy, those that resolve the regulator `η`
loop_dims(vac, η, tol::Real) = auto_loop_grid(vac, η, tol)
loop_dims(vac, η, grid::BZGrid) = grid.dims

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

# Accumulates the rank-one masses c_i y_i y_i† at bath energies xs[i], with y_i
# the columns of `Y`. Each mass is split linearly between its two nearest bins.
# The weight c carries the sign σ of the channel, and would carry its thermal
# factor at T > 0. Sorting the shares by bin makes each bin one rank-k update,
# of the scaled columns stored in the workspace `Z`, of size at least dim × 2n.
function accum_binned!(ρ::PairMeasure, xs, cs, Y; Z=zeros(ComplexF64, ρ.dim, 2length(xs)))
    σ = sign(sum(cs))
    all(c -> σ * c ≥ 0, cs) || error("Masses of one measure must share a sign")
    ts = xs / ρ.Δ
    js = floor.(Int, ts)
    bins = [js; js .+ 1]
    weights = abs.([cs .* (js .+ 1 - ts); cs .* (ts - js)])
    perm = sortperm(bins)
    Z = view(Z, :, eachindex(perm))
    Z .= view(Y, :, mod1.(perm, length(xs))) .* sqrt.(weights[perm])'
    bins = bins[perm]
    S = zeros(ComplexF64, ρ.dim, ρ.dim)
    for j in unique(bins)
        BLAS.herk!('U', 'N', σ, view(Z, :, searchsorted(bins, j)), 0.0, S)
        v = get!(() -> zeros(ComplexF64, packed_index(ρ.dim, ρ.dim)), ρ.bins, j)
        for n′ in 1:ρ.dim, n in 1:n′
            v[packed_index(n, n′)] += S[n, n′]
        end
    end
end

# Cauchy transform K(z) = Σ_j ρ_j / (z - jΔ), summed over the measures `ρs` and
# returned as a `dim×dim×length(zs)` array. Each Hermitian bin is transformed as
# the real and imaginary parts of its packed upper triangle, real series that
# are then unpacked into K. If `zs` is a range whose step is a whole number of
# bins, the sum is a discrete convolution, evaluated by FFT where the bins are
# dense enough to pay for it.
function cauchy_transform(ρs, zs)
    dim = first(ρs).dim
    K = zeros(ComplexF64, dim, dim, length(zs))
    for ρ in ρs
        isempty(ρ.bins) && continue
        y = if use_fft(ρ, zs)
            packed_transform_fft(ρ, zs, round(Int, real(ComplexF64(step(zs))) / ρ.Δ))
        else
            packed_transform(ρ, zs)
        end
        # The series S and A of each packed entry give K_ab += S + iA and
        # K_ba += S - iA
        for k in eachindex(zs), b in 1:dim, a in 1:b
            c = packed_index(a, b)
            (S, A) = (y[2c-1, k], y[2c, k])
            K[a, b, k] += S + im * A
            a == b || (K[b, a, k] += S - im * A)
        end
    end
    return K
end

# Whether the FFT pays for itself. The explicit sum costs nbins × nz per entry,
# the FFT (m + 2) N log N, at a similar cost per operation (measured on Apple M5).
function use_fft(ρ::PairMeasure, zs)
    zs isa AbstractRange && length(zs) > 1 || return false
    # A range of complex numbers may store its step in extended precision
    dz = ComplexF64(step(zs))
    iszero(imag(dz)) || return false
    m = round(Int, real(dz) / ρ.Δ)
    m > 0 && m * ρ.Δ ≈ real(dz) || return false
    (j0, j1) = extrema(keys(ρ.bins))
    N = length(zs) + cld(j1 - j0 + 1, m)
    return length(ρ.bins) * length(zs) > (m + 2) * N * log2(N)
end

# Transforms of the packed series of `ρ` at each of `zs`, as a matrix with one
# row per series, by explicit sum over bins.
function packed_transform(ρ::PairMeasure, zs)
    y = zeros(ComplexF64, 2packed_index(ρ.dim, ρ.dim), length(zs))
    for (j, v) in ρ.bins, (k, z) in enumerate(zs)
        c = 1 / (z - j * ρ.Δ)
        @views y[:, k] .+= c .* reinterpret(Float64, v)
    end
    return y
end

# As `packed_transform`, at zs[k+1] = z₀ + k m Δ for k = 0, …, n-1. Writing each
# bin index as j = j₀ + m i + r, with 0 ≤ r < m,
#
#     K_k = Σ_r Σ_i ρ_{j₀+mi+r} h^r_{k-i},   h^r_s = 1 / (z₀ - (j₀ + r)Δ + s m Δ),
#
# is a sum of m convolutions on the grid of `zs`. These are summed in Fourier
# space, so that a single inverse transform serves every r, and a circular
# convolution of length N ≥ n + ni - 1 is exact at k < n. The real series are
# convolved by real FFTs with the real and imaginary parts of h, a block of `nc`
# series at a time.
function packed_transform_fft(ρ::PairMeasure, zs, m; nc=64)
    (; Δ, bins) = ρ
    n = length(zs)
    ns = 2packed_index(ρ.dim, ρ.dim)
    (j0, j1) = extrema(keys(bins))
    ni = cld(j1 - j0 + 1, m)
    N = nextprod((2, 3, 5), n + ni - 1)

    # The series of each bin j0 + mi + r, at vs[mi + r + 1]
    empty = zeros(ComplexF64, ns ÷ 2)
    vs = [reinterpret(Float64, get(bins, j, empty)) for j in j0:j0+m*ni-1]

    h = [1 / (first(zs) - (j0 + r) * Δ + (k < n ? k : k - N) * m * Δ) for k in 0:N-1, r in 0:m-1]
    (ĥr, ĥi) = (FFTW.rfft(real(h), 1), FFTW.rfft(imag(h), 1))

    X = zeros(N, nc)
    X̂ = zeros(ComplexF64, N ÷ 2 + 1, nc)
    (Yr, Yi) = (similar(X̂), similar(X̂))
    (yr, yi) = (similar(X), similar(X))
    pf = FFTW.plan_rfft(X, 1)
    pb = FFTW.plan_brfft(X̂, N, 1)
    y = zeros(ComplexF64, ns, n)
    for cs in Iterators.partition(1:ns, nc)
        fill!(Yr, 0)
        fill!(Yi, 0)
        # Rows past ni stay zero, as the circular convolution requires
        for r in 0:m-1
            for i in 0:ni-1
                v = vs[m*i + r + 1]
                for (c′, c) in enumerate(cs)
                    X[i+1, c′] = v[c]
                end
            end
            mul!(X̂, pf, X)
            @. Yr += X̂ * $view(ĥr, :, r+1)
            @. Yi += X̂ * $view(ĥi, :, r+1)
        end
        mul!(yr, pb, Yr)
        mul!(yi, pb, Yi)
        for (c′, c) in enumerate(cs), k in 1:n
            y[c, k] = complex(yr[k, c′], yi[k, c′]) / N
        end
    end
    return y
end
