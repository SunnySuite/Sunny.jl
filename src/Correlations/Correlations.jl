#### BROADENING

abstract type AbstractBroadening end

struct Broadening{F <: Function} <: AbstractBroadening
    kernel :: F  # (ω_transfer - ω_excitation) -> intensity
end

struct NonstationaryBroadening{F <: Function} <: AbstractBroadening
    kernel :: F  # (ω_excitation, ω_transfer) -> intensity
end

function (b::Broadening)(ω1, ω2)
    b.kernel(ω2 - ω1)
end

function (b::NonstationaryBroadening)(ω1, ω2)
    b.kernel(ω1, ω2)
end

"""
    lorentzian(; fwhm)

Returns the function `(Γ/2) / (π*(x^2+(Γ/2)^2))` where `fwhm = Γ` is the full
width at half maximum.
"""
function lorentzian2(; fwhm)
    Γ = fwhm
    return Broadening(x -> (Γ/2) / (π*(x^2+(Γ/2)^2)))
end

"""
    gaussian(; {fwhm, σ})

Returns the function `exp(-x^2/2σ^2) / √(2π*σ^2)`. Either `fwhm` or `σ` must be
specified, where `fwhm = (2.355...) * σ` is the full width at half maximum.
"""
function gaussian2(; fwhm=nothing, σ=nothing)
    if sum(.!isnothing.((fwhm, σ))) != 1
        error("Either fwhm or σ must be specified.")
    end
    σ = Float64(@something σ (fwhm/2√(2log(2))))
    return Broadening(x -> exp(-x^2/2σ^2) / √(2π*σ^2))
end


#### Q-POINTS WITH OPTIONAL METADATA

abstract type AbstractQPoints end

struct QPoints <: AbstractQPoints
    qs :: Vector{Vec3}
end

struct QPath <: AbstractQPoints
    qs :: Vector{Vec3}
    xticks :: Tuple{Vector{Int64}, Vector{String}}
end

struct QGrid{N} <: AbstractQPoints
    qs :: Vector{Vec3}
    q0 :: Vec3
    Δqs :: NTuple{N, Vec3}
    lengths :: NTuple{N, Int}
end

function Base.convert(::Type{AbstractQPoints}, x::AbstractArray)
    return QPoints(collect(Vec3.(x)))
end

function Base.show(io::IO, qpts::AbstractQPoints)
    print(io, string(typeof(qpts)) * " ($(length(qpts.qs)) samples)")
end

function Base.show(io::IO, ::MIME"text/plain", qpts::QPath)
    printstyled(io, "QPath ($(length(qpts.qs)) samples)\n"; bold=true, color=:underline)
    println(io, "  " * join(qpts.xticks[2], " → "))
end


"""
    q_space_path(cryst::Crystal, qs, n; labels=nothing)

Returns a 1D path consisting of `n` wavevectors sampled piecewise-linearly
between the `qs`. Although the `qs` are provided in reciprocal lattice units
(RLU), consecutive samples are spaced uniformly in the global (inverse-length)
coordinate system. Optional `labels` can be associated with each special
q-point, and will be used in plotting functions.
"""
function q_space_path(cryst::Crystal, qs, n; labels=nothing)
    length(qs) >= 2 || error("Include at least two wavevectors in list qs.")
    qs = Vec3.(qs)
    # Displacement vectors in RLU
    dqs = qs[begin+1:end] - qs[begin:end-1]

    # Determine ms, the number of points in each segment. First point is placed
    # at the beginning of segment. Each m scales like dq in absolute units. The
    # total should be sum(ms) == n-1, anticipating a final point for qs[end].
    ws = [norm(cryst.recipvecs * dq) for dq in dqs]
    ms_ideal = (n - 1) .* ws / sum(ws)
    ms = round.(Int, ms_ideal)
    delta = sum(ms) - (n - 1)
    if delta < 0
        # add points where m < m_ideal
        idxs = sortperm(ms - ms_ideal; rev=false)[1:abs(delta)]
        ms[idxs] .+= 1
    elseif delta > 0
        # remove points where m > m_ideal
        idxs = sortperm(ms - ms_ideal; rev=true)[1:abs(delta)]
        ms[idxs] .-= 1
    end
    @assert sum(ms) == n - 1

    # Each segment should have at least one sample point
    any(iszero, ms) && error("Increase sample points n")

    # Linearly interpolate on each segment
    path = Vec3[]
    markers = Int[]
    for (i, m) in enumerate(ms)
        push!(markers, 1+length(path))
        for j in 0:m-1
            push!(path, qs[i] + (j/m)*dqs[i])
        end
    end
    push!(markers, 1+length(path))
    push!(path, qs[end])

    labels = @something labels fractional_vec3_to_string.(qs)
    xticks = (markers, labels)
    return QPath(path, xticks)
end


#### MEASUREMENT

# Op is the type of a local observable operator. Either a Vec3 (for :dipole
# mode, in which case the observable is `op⋅S`) or a HermitianC64 (for :SUN
# mode, in which case op is an N×N matrix).
struct Measurement{Op <: Union{Vec3, HermitianC64}, F, Ret}
    # TODO: Don't PACK?
    observables :: Array{Op, 5}          # (latsize, natoms, nobs)
    corr_pairs :: Vector{NTuple{2, Int}} # (ncorr)
    combiner :: F                        # (q::Vec3, obs) -> Ret

    # TODO: Default combiner will be SVector?
    function Measurement(observables::Array{Op, 5}, corr_pairs, combiner::F) where {Op, F}
        # Lift return type of combiner function to type-level
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        @assert isbitstype(Ret)
        return new{Op, F, Ret}(observables, corr_pairs, combiner)
    end
end

Base.eltype(::Measurement{Op, F, Ret}) where {Op, F, Ret} = Ret


function all_dipole_observables(sys::System{0}; apply_g)
    observables = zeros(Vec3, size(eachsite(sys))..., 3)
    for site in eachsite(sys)
        # Component α of observable is op⋅S = g[α,β] S[β]. Minus sign would
        # cancel because observables come in pairs.
        op = apply_g ? sys.gs[site] : Mat3(I)
        for α in 1:3
            observables[site, α] = op[α, :]
        end
    end
    return observables
end

function all_dipole_observables(sys::System{N}; apply_g) where {N}
    observables = Array{HermitianC64, 5}(undef, size(eachsite(sys))..., 3)
    for site in eachsite(sys)
        S = spin_matrices_of_dim(; N=sys.Ns[site])
        op = apply_g ? sys.gs[site]*S : S
        for α in 1:3
            observables[site, α] = op[α]
        end
    end
    return observables
end

"""
    DSSF(sys::System; apply_g=true)

Specify measurement of the dynamical spin structure factor (DSSF) or its
"instantaneous" variant. The "full" structure factor intensity
``\\mathcal{S}^{αβ}(𝐪,ω)`` is returned as a 3×3 matrix, with indices ``(α, β)``
denoting dipole components in Cartesian coordinates.

By default, the g-factor or tensor is applied at each site, yielding a
correlation between magnetic moments. Set `apply_g = false` to measure a true
spin-spin correlation.

Intended for use with [`intensities`](@ref), [`intensities_bands`](@ref), and
related functions. See also the Sunny documentation on [Structure Factor
Calculations](@ref) for more details.
"""
function DSSF(sys::System; apply_g=true)
    observables = all_dipole_observables(sys; apply_g)
    combiner(_, data) = SA[
        data[6]       data[5]       data[3]
        conj(data[5]) data[4]       data[2]
        conj(data[3]) conj(data[2]) data[1]
    ]
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    return Measurement(observables, corr_pairs, combiner)
end

"""
    DSSF_perp(sys::System; apply_g=true)

Specify measurement of the dynamical spin structure factor (DSSF). Like
[`DSSF`](@ref), but contracts the 3×3 structure factor matrix with
``(I-𝐪⊗𝐪/q²)``, which projects perpendicular to the direction of momentum
transfer ``𝐪``. The contracted structure factor can be interpreted as a
scattering intensity for an unpolarized neutron beam, up to constant scaling
factors. In the singular limit ``𝐪 → 0``, the contraction matrix is replaced by
its rotational average, ``(2/3) I``.

See also [`DSSF`](@ref).
"""
function DSSF_perp(sys::System; apply_g=true)
    observables = all_dipole_observables(sys; apply_g)
    function combiner(q, data)
        q2 = norm2(q)
        # Imaginary part cancels by symmetric contraction
        data = real.(SVector{6}(data))
        dssf = SA[
            data[6] data[5] data[3]
            data[5] data[4] data[2]
            data[3] data[2] data[1]
        ]
        tr_dssf = tr(dssf)
        # "S-perp" contraction matrix (1 - q⊗q/q²) appropriate to unpolarized
        # neutrons. In the limit q → 0, use (1 - q⊗q/q²) → 2/3, which
        # corresponds to a spherical average over uncorrelated data:
        # https://github.com/SunnySuite/Sunny.jl/pull/131
        return iszero(q2) ? (2/3)*tr_dssf : tr_dssf - dot(q, dssf, q) / q2
    end
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    return Measurement(observables, corr_pairs, combiner)
end

"""
    DSSF_trace(sys::System; apply_g=true)

Specify measurement of the dynamical spin structure factor (DSSF). Like
[`DSSF`](@ref), but returns only the trace of the 3×3 structure factor matrix.
This quantity can be useful for checking quantum sum rules.

See also [`DSSF`](@ref).
"""
function DSSF_trace(sys::System{N}; apply_g=true) where N
    observables = all_dipole_observables(sys; apply_g)
    combiner(_, data) = real(data[1] + data[2] + data[3])
    corr_pairs = [(3,3), (2,2), (1,1)]
    return Measurement(observables, corr_pairs, combiner)
end


#### INTENSITIES

abstract type AbstractIntensities end

struct BandIntensities{T} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: AbstractQPoints
    # Dispersion for each band
    disp :: Array{Float64, 2} # (nbands × nq)
    # Intensity data as Dirac-magnitudes
    data :: Array{T, 2} # (nbands × nq)
end

struct BroadenedIntensities{T} <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # Wavevectors in RLU
    qpts :: AbstractQPoints
    # Regular grid of energies
    energies :: Vector{Float64}
    # Convolved intensity data
    data :: Array{T, 2} # (nω × nq)
end

struct PowderIntensities <: AbstractIntensities
    # Original chemical cell
    crystal :: Crystal
    # q magnitudes in inverse length
    radii :: Vector{Float64}
    # Regular grid of energies
    energies :: Vector{Float64}
    # Convolved intensity data
    data :: Array{Float64, 2} # (nω × nq)
end

function Base.show(io::IO, res::AbstractIntensities)
    sizestr = join(size(res.data), "×")
    print(io, string(typeof(res)) * " ($sizestr elements)")
end


# Returns |1 + nB(ω)| where nB(ω) = 1 / (exp(βω) - 1) is the Bose function. See
# also `classical_to_quantum` which additionally "undoes" the classical
# Boltzmann distribution.
function thermal_prefactor(kT, ω)
    if iszero(kT)
        return ω >= 0 ? 1 : 0
    else
        @assert kT > 0
        return abs(1 / (1 - exp(-ω/kT)))
    end
end

function broaden!(data::AbstractMatrix{Ret}, bands::BandIntensities{Ret}, energies; kernel) where Ret
    energies = collect(energies)
    issorted(energies) || error("energies must be sorted")

    nω = length(energies)
    nq = size(bands.data, 2)
    (nω, nq) == size(data) || error("Argument data must have size ($nω×$nq)")

    cutoff = 1e-12 * Statistics.quantile(norm.(vec(bands.data)), 0.95)

    for iq in axes(bands.data, 2)
        for (ib, b) in enumerate(view(bands.disp, :, iq))
            norm(bands.data[ib, iq]) < cutoff && continue
            for (iω, ω) in enumerate(energies)
                data[iω, iq] += kernel(b, ω) * bands.data[ib, iq]
            end
            # If this broadening is a bottleneck, one can terminate when kernel
            # magnitude is small. This may, however, affect reference data used
            # in test suite.
            #=
                iω0 = searchsortedfirst(energies, b)
                for iω in iω0:lastindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
                for iω in iω0-1:-1:firstindex(energies)
                    ω = energies[iω]
                    x = kernel(b, ω) * bands.data[ib, iq]
                    data[iω, iq] += x
                    x < cutoff && break
                end
            =#
        end
    end

    return data
end

function broaden(bands::BandIntensities, energies; kernel)
    data = zeros(eltype(bands.data), length(energies), size(bands.data, 2))
    broaden!(data, bands, energies; kernel)
    return BroadenedIntensities(bands.crystal, bands.qpts, collect(energies), data)
end

function calculate_excitations!(V, H, swt::SpinWaveTheory, q)
    (; sys) = swt
    q_global = orig_crystal(sys).recipvecs * q
    q_reshaped = sys.crystal.recipvecs \ q_global

    if sys.mode == :SUN
        swt_hamiltonian_SUN!(H, swt, q_reshaped)
    else
        @assert sys.mode in (:dipole, :dipole_large_S)
        swt_hamiltonian_dipole!(H, swt, q_reshaped)
    end

    try
        return bogoliubov!(V, H)
    catch _
        error("Instability at wavevector q = $q")
    end
end

function excitations(swt::SpinWaveTheory, q)
    L = nbands(swt)
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    return calculate_excitations!(V, H, swt, q)
end

"""
    dispersion(swt::SpinWaveTheory, qpts)

Given a list of wavevectors `qpts` in reciprocal lattice units (RLU), returns
excitation energies for each band. The return value `ret` is 2D array, and
should be indexed as `ret[band_index, q_index]`.
"""
function dispersion(swt::SpinWaveTheory, qpts)
    qpts = convert(AbstractQPoints, qpts)
    return reduce(hcat, excitations.(Ref(swt), qpts.qs))
end

function localize_observable(v::Vec3, data::SWTDataDipole, site::Int)
    R = data.local_rotations[site]
    return R' * v
end

function localize_observable(A::HermitianC64, data::SWTDataSUN, site::Int)
    U = data.local_unitaries[:, :, site]
    return Hermitian(U' * A * U)
end

function localize_observables(obs::Array{Op, 5}, data) where {Op}
    Nsites = prod(size(obs)[1:4])
    Nobs = size(obs, 5)
    obs = reshape(obs, (Nsites, Nobs))
    ret = copy(obs)
    for α in 1:Nobs, site in 1:Nsites
        ret[site, α] = localize_observable(obs[site, α], data, site)
    end
    return ret
end


"""
    intensities_bands(swt::SpinWaveTheory, qpts; formfactors=nothing, measure)

Calculate spin wave excitation bands for a set of q-points in reciprocal space. TODO.
"""
function intensities_bands(swt::SpinWaveTheory, qpts; formfactors=nothing, measure::Measurement{Op, F, Ret}) where {Op, F, Ret}
    qpts = convert(AbstractQPoints, qpts)
    (; sys) = swt
    cryst = orig_crystal(sys)

    # Number of atoms in magnetic cell
    @assert sys.latsize == (1,1,1)
    Na = length(eachsite(sys))
    if Na != prod(size(measure.observables)[1:4])
        error("Size mismatch. Check that SpinWaveTheory and Measurement were built from same System.")
    end

    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Preallocation
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    Avec_pref = zeros(ComplexF64, Na)
    disp = zeros(Float64, L, Nq)
    intensity = zeros(Ret, L, Nq)

    # Temporary storage for pair correlations
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    Nobs = size(measure.observables, 5)
    obs_local_frame = localize_observables(measure.observables, swt.data)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sys.crystal)

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q
        disp[:, iq] .= calculate_excitations!(V, H, swt, q)

        for i in 1:Na
            r_global = global_position(sys, (1,1,1,i))
            Avec_pref[i] = exp(- im * dot(q_global, r_global))
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], norm2(q_global))
        end

        Avec = zeros(ComplexF64, Nobs)

        # Fill `intensity` array
        for band = 1:L
            fill!(Avec, 0)
            if sys.mode == :SUN
                N = sys.Ns[1]
                v = reshape(view(V, :, band), N-1, Na, 2)
                for i in 1:Na, μ in 1:Nobs
                    @views O = swt.data.observables_localized[:, :, μ, i]
                    @assert O ≈ obs_local_frame[i, μ]

                    for α in 1:N-1
                        Avec[μ] += Avec_pref[i] * (O[α, N] * v[α, i, 2] + O[N, α] * v[α, i, 1])
                    end
                end
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                (; sqrtS) = swt.data
                v = reshape(view(V, :, band), Na, 2)
                for i in 1:Na, μ in 1:Nobs
                    # @views O = swt.data.observables_localized[:, :, μ, i]
                    # @assert O ≈ - obs_local_frame[i, μ]'

                    # This is the Avec of the two transverse and one
                    # longitudinal directions in the local frame. (In the
                    # local frame, z is longitudinal, and we are computing
                    # the transverse part only, so the last entry is zero)
                    displacement_local_frame = SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                    Avec[μ] += Avec_pref[i] * (sqrtS[i]/sqrt(2)) * (obs_local_frame[i, μ]' * displacement_local_frame)[1]
                end
            end

            map!(corrbuf, measure.corr_pairs) do (α, β)
                Avec[α] * conj(Avec[β]) / Ncells
            end
            intensity[band, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return BandIntensities{Ret}(cryst, qpts, disp, intensity)
end

"""
    intensities(swt::SpinWaveTheory, qpts; energies, kernel, formfactors=nothing, measure)

Calculate spin wave intensities for a set of q-points in reciprocal space. TODO.
"""
function intensities(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, formfactors=nothing, measure::Measurement)
    return broaden(intensities_bands(swt, qpts; formfactors, measure), energies; kernel)
end


"""
    powder_average(f, cryst, radii, n; seed=0)

Calculate a powder-average over structure factor intensities. The `radii` have
units of inverse length, each defines a spherical shell in reciprocal space. The
[Fibonacci lattice](https://arxiv.org/abs/1607.04590) yields `n` points on the
shells, with quasi-uniformity. Sample points on different shells are
decorrelated through random rotations. A consistent random number `seed` will
yield reproducible results. The function `f` should accept a list of q-points
and call a variant of [`intensities`](@ref).

# Example
```julia
radii = range(0.0, 3.0, 200)
res = powder_average(cryst, radii, 500) do qs
    intensities(swt, qs; energies, kernel, measure)
end
plot_intensities(res)
```
"""
function powder_average(f, cryst, radii, n; seed=0)
    (; energies) = f([Vec3(0,0,0)])
    rng = Random.Xoshiro(seed)
    data = zeros(length(energies), length(radii))
    qs = reciprocal_space_shell(cryst, 1.0, n)
    
    for (i, radius) in enumerate(radii)
        R = Mat3(random_orthogonal(rng, 3)) * radius
        res = f(Ref(R) .* qs)
        data[:, i] = Statistics.mean(res.data; dims=2)
    end

    return PowderIntensities(cryst, radii, energies, data)
end



"""
    rotation_in_rlu(cryst::Crystal, (axis, angle))
    rotation_in_rlu(cryst::Crystal, R)

Returns a ``3×3`` matrix that rotates wavevectors in reciprocal lattice units
(RLU), with possible reflection. The input should be a representation of this
same rotation in global coordinates, i.e., a transformation of reciprocal-space
wavevectors in units of inverse length.
"""
function rotation_in_rlu end

function rotation_in_rlu(cryst::Crystal, (axis, angle))
    return rotation_in_rlu(cryst, axis_angle_to_matrix(axis, angle))
end

function rotation_in_rlu(cryst::Crystal, rotation::R) where {R <: AbstractMatrix}
    return inv(cryst.recipvecs) * Mat3(rotation) * cryst.recipvecs
end


"""
    domain_average(f, cryst, qpts; rotations, weights=nothing)

Calculate an average intensity for the reciprocal-space points `qpts` under a
discrete set of `rotations`. Rotations must be given in global Cartesian
coordinates, and will be converted via [`rotation_in_rlu`](@ref). Either
axis-angle or 3×3 rotation matrix representations can be used. An optional list
of `weights` allows for non-uniform weighting of each rotation. The function `f`
should accept a list of rotated q-points and call a variant of
[`intensities`](@ref).

# Example

```julia
# 0, 120, and 240 degree rotations about the global z-axis
rotations = [([0,0,1], n*(2π/3)) for n in 0:2]
weights = [1, 1, 1]
res = domain_average(cryst, qpts; rotations, weights) do qpts_rotated
    intensities(swt, qpts_rotated; energies, kernel, measure)
end
plot_intensities(res)
```
"""
function domain_average(f, cryst, qpts; rotations, weights=nothing)
    isempty(rotations) && error("Rotations must be nonempty list")
    weights = @something weights fill(1, length(rotations))
    length(rotations) == length(weights) || error("Rotations and weights must be same length")

    R0, Rs... = rotation_in_rlu.(Ref(cryst), rotations)
    w0, ws... = weights

    qpts = convert(AbstractQPoints, qpts)
    qs0 = copy(qpts.qs)

    qpts.qs .= Ref(R0) .* qs0
    res = f(qpts)
    res.data .*= w0

    for (R, w) in zip(Rs, ws)
        qpts.qs .= Ref(R) .* qs0
        res.data .+= w .* f(qpts).data
    end

    qpts.qs .= qs0
    res.data ./= sum(weights)
    return res
end
