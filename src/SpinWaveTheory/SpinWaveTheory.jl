###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################

struct SWTDataDipole
    R_mat  :: Vector{Mat3}             # SO(3) rotation to align the quantization axis
    c_coef :: Vector{StevensExpansion} # Stevens operator coefficents
end

struct SWTDataSUN
    s̃_mat :: Array{ComplexF64, 4}  # Dipole operators
    T̃_mat :: Array{ComplexF64, 3}  # Single-ion anisos
    Q̃_mat :: Array{ComplexF64, 4}  # Quadrupolar operators
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
    cellsize_mag = cell_dimensions(sys) * diagm(collect(sys.latsize))
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
    Nₘ = natoms(sys.crystal)

    s_mat_N = spin_matrices(; N)

    # we support the biquad interactions now in the :dipole mode
    # we choose a particular basis of the nematic operators listed in Appendix B of *Phys. Rev. B 104, 104409*
    Q_mat_N = Vector{Matrix{ComplexF64}}(undef, 5)
    Q_mat_N[1] = -(s_mat_N[1] * s_mat_N[3] + s_mat_N[3] * s_mat_N[1])
    Q_mat_N[2] = -(s_mat_N[2] * s_mat_N[3] + s_mat_N[3] * s_mat_N[2])
    Q_mat_N[3] = s_mat_N[1] * s_mat_N[1] - s_mat_N[2] * s_mat_N[2]
    Q_mat_N[4] = s_mat_N[1] * s_mat_N[2] + s_mat_N[2] * s_mat_N[1]
    Q_mat_N[5] = √3 * s_mat_N[3] * s_mat_N[3] - 1/√3 * S * (S+1) * I

    U_mat = Matrix{ComplexF64}(undef, N, N)

    s̃_mat = Array{ComplexF64, 4}(undef, N, N, 3, Nₘ)
    T̃_mat = Array{ComplexF64, 3}(undef, N, N, Nₘ)
    Q̃_mat = zeros(ComplexF64, N, N, 5, Nₘ)

    for atom = 1:Nₘ
        U_mat[:, 1] = sys.coherents[1, 1, 1, atom]
        U_mat[:, 2:N] = nullspace(U_mat[:, 1]')
        for μ = 1:3
            s̃_mat[:, :, μ, atom] = Hermitian(U_mat' * s_mat_N[μ] * U_mat)
        end
        for ν = 1:5
            Q̃_mat[:, :, ν, atom] = Hermitian(U_mat' * Q_mat_N[ν] * U_mat)
        end
        T̃_mat[:, :, atom] = Hermitian(U_mat' * sys.interactions_union[atom].onsite.matrep * U_mat)
    end

    return SWTDataSUN(s̃_mat, T̃_mat, Q̃_mat)
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
        c = rotate_operator(sys.interactions_union[atom].onsite.stvexp, R)

        push!(Rs, R)
        push!(cs, c)
    end

    return SWTDataDipole(Rs, cs)
end
