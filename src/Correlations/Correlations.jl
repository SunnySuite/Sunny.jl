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

function lorentzian2(; fwhm)
    Γ = fwhm
    return Broadening(x -> (Γ/2) / (π*(x^2+(Γ/2)^2)))
end


#### Q-POINTS WITH OPTIONAL METADATA

abstract type AbstractQPoints end

struct QPoints <: AbstractQPoints
    qs :: Vector{Vec3}
end

struct QPath <: AbstractQPoints
    qs :: Vector{Vec3}
    xticks :: Vector{Tuple{Int, String}}
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


"""
    q_space_path(cryst::Crystal, qs, Δq)

Returns a 1D path in q-space that samples linearly between the provided
wavevectors `qs`, with approximate displacements Δq between steps. Both
`qs` and `Δq` must be provided in reciprocal lattice units (RLU) for the
given crystal.
"""
function q_space_path(cryst::Crystal, qs, Δq)
    @assert length(qs) >= 2 "The list `qs` should include at least two wavevectors."
    qs = Vec3.(qs)

    Δq_absolute = minimum(norm.(eachcol(cryst.recipvecs))) * Δq

    path = Vec3[]
    markers = Int[]
    for i in 1:length(qs)-1
        push!(markers, length(path)+1)
        q1, q2 = qs[i], qs[i+1]
        dist = norm(cryst.recipvecs * (q1 - q2))
        npoints = ceil(Int, dist/Δq_absolute)
        for n in 1:npoints
            push!(path, (1 - (n-1)/npoints)*q1 + (n-1)*q2/npoints)
        end
    end
    push!(markers, length(path)+1)
    push!(path, qs[end])

    labels = map(qs) do q
        # "[" * join(number_to_math_string.(q), ",") * "]"
        fractional_vec3_to_string(q)
    end
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
        # Component α of observable is op⋅S = -g[α,β] S[β]
        M = apply_g ? sys.gs[site] : Mat3(I)  # TODO: Sign
        for α in 1:3
            observables[site, α] = M[α, :]
        end
    end
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    return observables, corr_pairs
end

function all_dipole_observables(sys::System{N}; apply_g) where {N}
    observables = Array{HermitianC64, 5}(undef, size(eachsite(sys))..., 3)
    for site in eachsite(sys)
        S = spin_matrices_of_dim(; N=sys.Ns[site])
        M = apply_g ? sys.gs[site]*S : S  # TODO: Sign
        for α in 1:3
            observables[site, α] = M[α]
        end
    end
    corr_pairs = [(3,3), (2,3), (1,3), (2,2), (1,2), (1,1)]
    return observables, corr_pairs
end


function DSSF(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    combiner(_, data) = SA[
        data[6]       data[5]       data[3]
        conj(data[5]) data[4]       data[2]
        conj(data[3]) conj(data[2]) data[1]
    ]
    return Measurement(observables, corr_pairs, combiner)
end

function DSSF_perp(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    function combiner(q, data)
        q2 = norm2(q)
        # Imaginary part cancels by symmetric contraction
        data = real.(SVector{6}(data))
        dssf = SA[
            data[6] data[5] data[3]
            data[5] data[4] data[2]
            data[3] data[2] data[1]
        ]
        return tr(dssf) - (q' * dssf * q) / (q2 + 1e-14)
        # TODO: Check how SpinW regularizes, and consider also "Mourigal limit",
        # https://github.com/SunnySuite/Sunny.jl/pull/131
    end
    return Measurement(observables, corr_pairs, combiner)
end

function DSSF_trace(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    combiner(_, data) = real(data[1] + data[4] + data[6])
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

# TODO: "excitations!" ?
function calculate_quasiparticles!(V, H, swt::SpinWaveTheory, q)
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

# TODO: excitations()
function dispersion2(swt::SpinWaveTheory, q)
    L = nbands(swt)
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    return calculate_quasiparticles!(V, H, swt, q)
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

# TODO: measure=nothing
function intensities2(swt::SpinWaveTheory, qpts; formfactors=nothing, measure::Measurement{Op, F, Ret}) where {Op, F, Ret}
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
        disp[:, iq] .= calculate_quasiparticles!(V, H, swt, q)

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
                v = reshape(view(V, :, band), Na, 2)
                for i in 1:Na
                    sqrtS = sqrt(sys.κs[i])
                    for μ in 1:Nobs
                        @views O = swt.data.observables_localized[:, :, μ, i]
                        # @assert O ≈ - obs_local_frame[i, μ]'

                        # This is the Avec of the two transverse and one
                        # longitudinal directions in the local frame. (In the
                        # local frame, z is longitudinal, and we are computing
                        # the transverse part only, so the last entry is zero)
                        displacement_local_frame = SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                        Avec[μ] += Avec_pref[i] * (sqrtS/sqrt(2)) * (obs_local_frame[i, μ]' * displacement_local_frame)[1]
                    end
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


function intensities_broadened2(swt::SpinWaveTheory, qpts, energies; kernel::B, formfactors=nothing, measure::Measurement) where {B <: AbstractBroadening}
    all(energies .>= 0) || error("Energies must be non-negative")

    qpts = convert(AbstractQPoints, qpts)
    energies = collect(energies)
    (; qs) = qpts
    bands = intensities2(swt, qpts; formfactors, measure)

    nω = length(energies)
    nq = length(qs)
    data = zeros(eltype(measure), nω, nq)

    for iq in eachindex(qs)
        for (ib, b) in enumerate(view(bands.disp, :, iq))
            for (iω, ω) in enumerate(energies)
                data[iω, iq] += kernel(b, ω) * bands.data[ib, iq]
            end
        end
    end

    return BroadenedIntensities(bands.crystal, qpts, energies, data)
end

