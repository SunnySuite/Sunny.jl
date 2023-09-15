###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################

struct SWTDataDipole
    R_mat  :: Vector{Mat3}             # SO(3) rotation to align the quantization axis
    c_coef :: Vector{StevensExpansion} # Stevens operator coefficents
end

struct SWTDataSUN
    dipole_operators :: Array{ComplexF64, 4}
    onsite_operator :: Array{ComplexF64, 3}  # Single-ion anisotropy
    quadrupole_operators :: Array{ComplexF64, 4}
end

"""
    SpinWaveTheory(sys, energy_ϵ::Float64=1e-8, energy_tol=1e-6)

Constructs an object to perform linear spin wave theory. Use it with
[`dispersion`](@ref) and [`dssf`](@ref) functions.

The optional parameter `energy_ϵ` adds a small positive shift to the diagonal of
the dynamical matrix ``D`` to avoid numerical issues with zero-energy
quasi-particle modes. The optional parameter `energy_tol` relaxes the check on
the imaginary part of the eigenvalues.
"""
struct SpinWaveTheory
    sys        :: System
    data       :: Union{SWTDataDipole, SWTDataSUN}
    energy_ϵ   :: Float64          # Energy shift applied to dynamical matrix prior to Bogoliubov transformation
    energy_tol :: Float64          # Energy tolerance for maximal imaginary part of quasiparticle energies

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    # dipole_corrs :: Bool                                  # Whether using all correlations from dipoles 
    # observables  :: Array{ComplexF64, 3}                  # Operators corresponding to observables
    # idxinfo      :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)
end

function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, energy_tol::Float64=1e-6) where N
    if !isnothing(sys.ewald)
        error("SpinWaveTheory does not yet support long-range dipole-dipole interactions.")
    end

    # Reshape into single unit cell. A clone will always be performed, even if
    # no reshaping happens.
    cellsize_mag = cell_shape(sys) * diagm(collect(sys.latsize))
    sys = reshape_supercell_aux(sys, (1,1,1), cellsize_mag)

    # Rotate local operators to quantization axis
    data = sys.mode == :SUN ? swt_data_sun(sys) : swt_data_dipole(sys)

    return SpinWaveTheory(sys, data, energy_ϵ, energy_tol)
end


function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    # modename = swt.dipole_corrs ? "Dipole correlations" : "Custom correlations"
    modename = "Dipole correlations"
    printstyled(io, "SpinWaveTheory [$modename]\n"; bold=true, color=:underline)
    println(io, "Atoms in magnetic supercell: $(natoms(swt.sys.crystal))")
end

function num_bands(swt::SpinWaveTheory)
    (; sys) = swt
    nflavors = sys.mode == :SUN ? sys.Ns[1]-1 : 1
    return nflavors * natoms(sys.crystal)
end


# Convert 3-vector from the Cartesian frame to the spherical frame
function dipole_to_angles(dipoles)
    r = norm(dipoles)
    @assert r > 1e-7
    θ = acos(dipoles[3] / r)
    ϕ = atan(dipoles[2], dipoles[1])
    return θ, ϕ
end

# Given q in reciprocal lattice units (RLU) for the original crystal, return a
# q_reshaped in RLU for the possibly-reshaped crystal.
function to_reshaped_rlu(sys::System{N}, q) where N
    return sys.crystal.recipvecs \ (orig_crystal(sys).recipvecs * q)
end

# Compute SU(N) generators in the local reference frame (for :SUN mode). DD:
# Redo this using existing operator rotation facilities.
function swt_data_sun(sys::System{N}) where N
    S = (N-1)/2
    n_magnetic_atoms = natoms(sys.crystal)

    dipole_operators = spin_matrices(; N)
	Sx, Sy, Sz = dipole_operators

    # we support the biquad interactions now in the :dipole mode
    # we choose a particular basis of the nematic operators listed in Appendix B of *Phys. Rev. B 104, 104409*
    quadrupole_operators = Vector{Matrix{ComplexF64}}(undef, 5)
    quadrupole_operators[1] = -(Sx * Sz + Sz * Sx)
    quadrupole_operators[2] = -(Sy * Sz + Sz * Sy)
    quadrupole_operators[3] = Sx * Sx - Sy * Sy
    quadrupole_operators[4] = Sx * Sy + Sy * Sx
    quadrupole_operators[5] = √3 * Sz * Sz - 1/√3 * S * (S+1) * I

    local_quantization_basis = Matrix{ComplexF64}(undef, N, N)

    dipole_operators_localized = Array{ComplexF64, 4}(undef, N, N, 3, n_magnetic_atoms)
    onsite_operator_localized = Array{ComplexF64, 3}(undef, N, N, n_magnetic_atoms)
    quadrupole_operators_localized = zeros(ComplexF64, N, N, 5, n_magnetic_atoms)

    for atom = 1:n_magnetic_atoms
	    # First axis of local quantization basis is along the 
		# ground-state polarization axis
        local_quantization_basis[:, 1] = sys.coherents[1, 1, 1, atom]
		
		# Remaining axes are arbitrary but mutually orthogonal
		# and orthogonal to the first axis
        local_quantization_basis[:, 2:N] = nullspace(local_quantization_basis[:, 1]')
		
        for μ = 1:3
            dipole_operators_localized[:, :, μ, atom] = Hermitian(local_quantization_basis' * dipole_operators[μ] * local_quantization_basis)
        end
        for ν = 1:5
            quadrupole_operators_localized[:, :, ν, atom] = Hermitian(local_quantization_basis' * quadrupole_operators[ν] * local_quantization_basis)
        end
        onsite_operator_localized[:, :, atom] = Hermitian(local_quantization_basis' * sys.interactions_union[atom].onsite * local_quantization_basis)
    end

    return SWTDataSUN(dipole_operators_localized, onsite_operator_localized, quadrupole_operators_localized)
end

# Compute Stevens coefficients in the local reference frame
function swt_data_dipole(sys::System{0})
    cs = StevensExpansion[]
    Rs = Mat3[]

    for atom in 1:natoms(sys.crystal)
        # SO(3) rotation that aligns the quantization axis. Important: since we
        # will project out bosons that correspond to multipolar fluctuations,
        # therefore we use the explicit matrix to get rid of any ambiguity.
        #
        # As a unitary, U = exp(-i ϕ Sz) exp(-i θ Sy)
        θ, ϕ = dipole_to_angles(sys.dipoles[1,1,1,atom])
        R = SA[-sin(ϕ) -cos(ϕ)*cos(θ) cos(ϕ)*sin(θ);
                cos(ϕ) -sin(ϕ)*cos(θ) sin(ϕ)*sin(θ);
                0.0     sin(θ)        cos(θ)]
        # Rotated Stevens expansion
        c = rotate_operator(sys.interactions_union[atom].onsite, R)

        push!(Rs, R)
        push!(cs, c)
    end

    return SWTDataDipole(Rs, cs)
end
