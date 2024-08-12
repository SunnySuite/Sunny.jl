#### Q-POINTS

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

"""
    blume_maleev(q_i, q)

Given an incident momentum `q_i` and a momentum transfer vector to the sample
`q`, this function returns three orthogonal basis vectors `(ex, ey, ez)` that
define the Blume-Maleev polarization axis system. All input and output vectors
are given in the global Cartesian system for reciprocal space (inverse length
units).

The three Blume-Maleev basis vectors are constructed as follows:

```julia
ex = normalize(q)       # parallel to q
ez = normalize(q × q_i) # normal to the scattering plane
ey = normalize(ez × q)  # perpendicular to q, in the scattering plane
```
"""
function blume_maleev(q_i, q)
    q_i, q = Vec3.((q_i, q))
    ex = normalize(q)
    ez = normalize(q × q_i)
    ey = normalize(ez × q)
    return (ex, ey, ez)
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


function broaden!(data::AbstractMatrix{Ret}, bands::BandIntensities{Ret}; energies, kernel) where Ret
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

function broaden(bands::BandIntensities; energies, kernel)
    data = zeros(eltype(bands.data), length(energies), size(bands.data, 2))
    broaden!(data, bands; energies, kernel)
    return BroadenedIntensities(bands.crystal, bands.qpts, collect(energies), data)
end


#### ROTATIONAL AVERAGING

"""
    powder_average(f, cryst, radii, n; seed=0)

Calculate a powder-average over structure factor intensities. The `radii`, with
units of inverse length, define spherical shells in reciprocal space. The
[Fibonacci lattice](https://arxiv.org/abs/1607.04590) yields `n` points on the
sphere, with quasi-uniformity. Sample points on different shells are
decorrelated through random rotations. A consistent random number `seed` will
yield reproducible results. The function `f` should accept a list of q-points
and call a variant of [`intensities`](@ref).

# Example
```julia
radii = range(0.0, 3.0, 200)
res = powder_average(cryst, radii, 500) do qs
    intensities(swt, qs; energies, kernel)
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
res = domain_average(cryst, path; rotations, weights) do path_rotated
    intensities(swt, path_rotated; energies, kernel)
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
