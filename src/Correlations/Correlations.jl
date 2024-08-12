# abstract type Broadening end

# struct OneParamBroadening
# end

# struct EnergyDependentBroadening
# end


# * Op: The type of a local observable operator. Either a Vec3 (for :dipole
#   mode, in which case the observable is `op⋅S`) or a HermitianC64 (for :SUN
#   mode, in which case op is already an N×N matrix).
# * Scalar: Type of scalar expectation value ⟨Op⟩. Either ComplexF64 or Float64.
#   The latter retains only the real part.
# * Ret: Type of the return value after running `contractor`. Must be `isbits`,
#   but many scalars can be packed into the single return value.
struct Measurement{Op, Scalar, Ret}
    observables :: Array{Op, 5} # (na, latsize, nobs)
    pair_indices :: Vector{NTuple{2, Int}} # (ncorr)
    # A function that takes (q::Vec3, data::Vector{Scalar}) and produces a
    # contracted result Ret
    contractor :: Function
end

function DSSF_observables(sys::System{N}; apply_g) where N

end

function DSSF_observables(sys::System{0}; apply_g)

end

function DSSF(sys::System{N}; apply_g=true)
    observables, pair_indices = DSSF_observables(sys; apply_g)
    contractor(_, data) = SA[
        data[1]       data[2]       data[4]
        conj(data[2]) data[3]       data[5]
        conj(data[4]) conj(data[5]) data[6]
    ]
    Ret = Hermitian{ComplexF64, SMatrix{3, 3, ComplexF64, 9}}
    return Measurement{eltype(observables), ComplexF64, Ret}(observables, pair_indices, contractor)
end

function DSSF_trace(sys::System{N}; apply_g=true)
    observables, pair_indices = DSSF_observables(sys; apply_g)
    contractor(_, data) = data[1] + data[3] + data[6]
    return Measurement{eltype(observables), Float64, Float64}(observables, pair_indices, contractor)
end

function DSSF_perp(sys::System{N}; apply_g=true)
    observables, pair_indices = DSSF_observables(sys; apply_g)
    function contractor(q, data)
        data[1] + data[3] + data[6]
    end
    return Measurement{eltype(observables), Float64, Float64}(observables, pair_indices, contractor)
end


abstract type AbstractIntensities end

struct BandIntensities{T} <: AbstractIntensities
    # Dispersion for each band
    disp :: Array{Float64, 2} # (nbands × nq)
    # Wavevectors in absolute units
    qs_abs :: Vector{Vec3}
    # Intensity data as Dirac-magnitudes
    data :: Array{T, 3} # (ncorr × nbands × nq)
end

struct BroadenedIntensities{T} <: AbstractIntensities
    # Regular grid of energies
    energies :: AbstractRange
    # Wavevectors in absolute units
    qs_abs :: Vector{Vec3}
    # Integrated intensity over Δω
    data :: Array{T, 3} # (ncorr × nω × nq)
end

# An function-like object that can be applied to momentum and intensity data (q,
# dat), and carries its return type statically.
struct IntensityContractor{T, F}
    f::F
end
(c::IntensityContractor)(q, dat) = c.f(q, dat)


const DssfFull = begin
    f = function(_, dat)
        @assert length(dat) == 6
        Hermitian(SA[dat[1]       dat[2]       dat[4]
                     conj(dat[2]) dat[3]       dat[5]
                     conj(dat[4]) conj(dat[5]) dat[6]])
    end
    HMat3 = Hermitian{ComplexF64, SMatrix{3, 3, ComplexF64, 9}}
    IntensityContractor{HMat3, typeof(f)}(f)
end

const DssfTrace = begin
    f = (q, dat) -> tr(DssfFull(q, dat))
    IntensityContractor{Float64, typeof(f)}(f)
end

const DssfPerp = begin
    f = function(q, dat)
        if iszero(q)
            DssfTrace(q, dat)
        else
            A = real(DssfFull(q, dat))
            B = Hermitian(I - q * q' / norm2(q))
            dot(A, B) # == tr(A' * B)
        end
    end
    IntensityContractor{Float64, typeof(f)}(f)
end



function intensities(swt::SpinWaveTheory, qs, formfactors; contractor)
    contractor = let
        if isnothing(contractor) || typeof(contractor) == IntensityContractor
            contractor
        elseif contractor == :full
            DssfFull
        elseif contractor == :trace
            DssfTrace
        elseif contractor == :perp
            DssfPerp
        else
            error("Unknown contractor")
        end
    end
    intensities_aux(swt, qs, formfactors, contractor)
end

function intensities_aux(swt::SpinWaveTheory, qs, formfactors, contractor::Union{Nothing, Contractor{T, F}}) where {T, F}

    # Number of atoms in magnetic cell
    @assert sys.latsize = (1,1,1)
    Na = natoms(sys.crystal)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(orig_crystal(sys))
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qs)
    Ncorr = length(observables.correlations)

    # Preallocation
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)
    Avec_pref = zeros(ComplexF64, Na)
    disp = zeros(Float64, L, Nq)
    intensity = isnothing(contractor) ? zeros(ComplexF64, Ncorr, L, Nq) : zeros(T, L, Nq)
    corrbuf = zeros(ComplexF64, Ncorr)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, swt.sys.crystal)
    
    for (iq, q) in enumerate(qs)
        q_absolute = orig_crystal(swt.sys).recipvecs * q
        q_reshaped = swt.sys.crystal.recipvecs \ q

        if sys.mode == :SUN
            swt_hamiltonian_SUN!(H, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H, swt, q_reshaped)
        end

        try
            disp[:, iq] .= bogoliubov!(V, H)
        catch _
            error("Instability at wavevector q = $q")
        end

        for i in 1:Na
            Avec_pref[i] = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)
        end

        Avec = zeros(ComplexF64, num_observables(observables))
        
        # Fill `intensity` array
        for band = 1:L
            fill!(Avec, 0)
            if sys.mode == :SUN
                (; observables_localized) = data::SWTDataSUN
                N = sys.Ns[1]
                v = reshape(view(V, :, band), N-1, Na, 2)
                for i in 1:Na
                    for μ in 1:num_observables(observables)
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

            if isnothing(contractor)
                intensity[:, band, iq] .= corrbuf
            else
                intensity[band, iq] = contractor(q, corrbuf)
            end
        end
    end

    BandIntensities{Float64}(
        disp, # (nbands × nq)
        [orig_crystal(swt.sys).recipvecs * q for q in qs], # qs_abs
        intensity, # (ncorr × nbands × nq)
    )
end
