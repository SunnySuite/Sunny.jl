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
    # A function that takes (q::Vec3, elems) and produces a combined result
    combiner :: F

    function Measurement(observables::Array{Op, 5}, corr_pairs, combiner::F) where {Op, F}
        # Lift return type of combiner function to type-level
        Ret = only(Base.return_types(combiner, (Vec3, Vector{ComplexF64})))
        @assert isbitstype(Ret)
        return new{Op, F, Ret}(observables, corr_pairs, combiner)
    end
end

function all_dipole_observables(sys::System; apply_g)
    observables = zeros(Vec3, size(eachsite(sys))..., 3)
    for site in eachsite(sys)
        M = apply_g ? -sys.gs[site] : Mat3(I)
        for α in 1:3
            observables[site, α] = M[α, :]
        end
    end
    corr_pairs = [(1,1), (1,2), (2,2), (1,3), (2,3), (3,3)]
    return observables, corr_pairs
end

function DSSF(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    combiner(_, data) = SA[
        data[1]       data[2]       data[4]
        conj(data[2]) data[3]       data[5]
        conj(data[4]) conj(data[5]) data[6]
    ]
    return Measurement(observables, corr_pairs, combiner)
end

function DSSF_trace(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    data = real.(SVector{6}(data))
    combiner(_, data) = data[1] + data[3] + data[6]
    return Measurement(observables, corr_pairs, combiner)
end

function DSSF_perp(sys::System{N}; apply_g=true) where N
    observables, corr_pairs = all_dipole_observables(sys; apply_g)
    function combiner(q, data)
        q2 = norm2(q)
        data = real.(SVector{6}(data))
        return data[1] + data[3] + data[6] +
            - (data[1]*q[1]^2 + data[3]*q[2]^2 + data[6]*q[3]^2) / q2 + 
            - 2 * (data[2]*q[1]*q[2] + data[4]*q[1]*q[3] + data[5]*q[2]*q[3]) / q2
        # A = [data[1] data[2] data[4]
        #      data[2] data[3] data[5]
        #      data[4] data[5] data[6]]
        # B = Hermitian(I - q * q' / norm2(q))
        # return tr(A' * B)
    end
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
    (; data, observables) = swt

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

        nobs = size(measure.observables, 5)
        Avec = zeros(ComplexF64, nobs)

        # Fill `intensity` array
        for band = 1:L
            fill!(Avec, 0)
            if sys.mode == :SUN
                (; observables_localized) = data::SWTDataSUN
                N = sys.Ns[1]
                v = reshape(view(V, :, band), N-1, Na, 2)
                for i in 1:Na
                    for μ in 1:nobs
                        @views O = observables_localized[:, :, μ, i]
                        for α in 1:N-1
                            Avec[μ] += Avec_pref[i] * (O[α, N] * v[α, i, 2] + O[N, α] * v[α, i, 1])
                        end
                    end
                end
            else
                (; observables_localized, sqrtS) = data::SWTDataDipole
                @assert sys.mode in (:dipole, :dipole_large_S)
                v = reshape(view(V, :, band), Na, 2)
                for i in 1:Na
                    for μ in 1:num_observables(observables)
                        # This is the Avec of the two transverse and one
                        # longitudinal directions in the local frame. (In the
                        # local frame, z is longitudinal, and we are computing
                        # the transverse part only, so the last entry is zero)
                        displacement_local_frame = SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]

                        # Note that O_local_frame has already been right
                        # multiplied by data.local_rotations[i]
                        @views O_local_frame = observables_localized[:,:,μ,i]
                        Avec[μ] += Avec_pref[i] * (sqrtS[i]/sqrt(2)) * (O_local_frame * displacement_local_frame)[1]
                    end
                end
            end

            for (c, ic) in observables.correlations
                (α, β) = c.I
                corrbuf[ic] = Avec[α] * conj(Avec[β]) / Ncells
            end
            intensity[band, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return BandIntensities{Ret}(qs_global, disp, intensity)
end
