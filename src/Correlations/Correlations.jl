# abstract type Broadening end

# struct OneParamBroadening
# end

# struct EnergyDependentBroadening
# end

# Op is the type of a local observable operator. Either a Vec3 (for :dipole
# mode, in which case the observable is `op⋅S`) or a HermitianC64 (for :SUN
# mode, in which case op is an N×N matrix).
struct Measurement{Op, F, Ret}
    observables :: Array{Op, 5}          # (latsize, natoms, nobs)
    corr_pairs :: Vector{NTuple{2, Int}} # (ncorr)
    combiner :: F                        # (q::Vec3, obs) -> Ret

    function Measurement(observables::Array{Op, 5}, corr_pairs, combiner::F) where {Op, F}
        # Lift return type of combiner function to type-level
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        @assert isbitstype(Ret)
        return new{Op, F, Ret}(observables, corr_pairs, combiner)
    end
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


function all_dipole_observables(sys::System{0}; apply_g)
    observables = zeros(Vec3, size(eachsite(sys))..., 3)
    for site in eachsite(sys)
        # Component α of observable is op⋅S = -g[α,β] S[β]
        M = apply_g ? -sys.gs[site] : Mat3(I)
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
        M = apply_g ? -sys.gs[site]*S : S
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
        return tr(dssf) - (q' * dssf * q) / q2
    end
    return Measurement(observables, corr_pairs, combiner)
end

function DSSF_trace(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    combiner(_, data) = real(data[1] + data[4] + data[6])
    return Measurement(observables, corr_pairs, combiner)
end


abstract type AbstractIntensities end

struct BandIntensities{T} <: AbstractIntensities
    # Wavevectors in global (inverse length) units
    qs_global :: Vector{Vec3}
    # Dispersion for each band
    disp :: Array{Float64, 2} # (nbands × nq)
    # Intensity data as Dirac-magnitudes
    data :: Array{T, 2} # (nbands × nq)
end

struct BroadenedIntensities{T} <: AbstractIntensities
    # Wavevectors in global (inverse length) units
    qs_global :: Vector{Vec3}
    # Regular grid of energies
    energies :: AbstractRange
    # Integrated intensity over Δω
    data :: Array{T, 2} # (nω × nq)
end


# TODO: Move measure into SWT
function intensities2(swt::SpinWaveTheory, qs; formfactors=nothing, measure::Measurement{Op, F, Ret}) where {Op, F, Ret}
    (; sys) = swt

    # Number of atoms in magnetic cell
    @assert sys.latsize == (1,1,1)
    Na = natoms(sys.crystal)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(orig_crystal(sys))
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qs)
    # Wavevectors in global (inverse length) units
    qs_global = [orig_crystal(sys).recipvecs * q for q in qs]

    # Preallocation
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)
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
    
    for iq in 1:Nq
        q_global = qs_global[iq]
        q_reshaped = sys.crystal.recipvecs \ q_global

        if sys.mode == :SUN
            swt_hamiltonian_SUN!(H, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H, swt, q_reshaped)
        end

        try
            disp[:, iq] .= bogoliubov!(V, H)
        catch _
            error("Instability at wavevector q = $(qs[iq])")
        end

        for i in 1:Na
            Avec_pref[i] = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))
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
                    @assert O ≈ - obs_local_frame[i, μ]

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
                        @assert O ≈ - obs_local_frame[i, μ]'

                        # This is the Avec of the two transverse and one
                        # longitudinal directions in the local frame. (In the
                        # local frame, z is longitudinal, and we are computing
                        # the transverse part only, so the last entry is zero)
                        displacement_local_frame = SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                        Avec[μ] += Avec_pref[i] * (sqrtS/sqrt(2)) * (obs_local_frame[i, μ]' * displacement_local_frame)[1]
                    end
                end
            end

            for (i, (α, β)) in enumerate(measure.corr_pairs)
                corrbuf[i] = Avec[α] * conj(Avec[β]) / Ncells
            end
            intensity[band, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return BandIntensities{Ret}(qs_global, disp, intensity)
end


function intensities_broadened2(swt::SpinWaveTheory, qs, energies::StepRangeLen; formfactors=nothing, measure::Measurement{Op, F, Ret}) where {Op, F, Ret}

end