###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################

"""
    SpinWaveTheory(sys, energy_ϵ::Float64=1e-8, energy_tol=1e-6)

**Experimental**. Constructs an object to perform linear spin wave theory. Use
it with [`dispersion`](@ref) and [`dssf`](@ref) functions.

The optional parameter `energy_ϵ` adds a small positive shift to the diagonal of
the dynamical matrix ``D`` to avoid numerical issues with zero-energy
quasi-particle modes. The optional parameter `energy_tol` relaxes the check on
the imaginary part of the eigenvalues.
"""
struct SpinWaveTheory
    sys   :: System
    s̃_mat :: Array{ComplexF64, 4}  # Dipole operators
    T̃_mat :: Array{ComplexF64, 3}  # Single-ion anisos
    Q̃_mat :: Array{ComplexF64, 4}  # Quadrupolar operators
    c′_coef :: Vector{StevensExpansion} # Stevens operator coefficents (for dipole mode)
    R_mat   :: Vector{Mat3}        # SO(3) rotation to align the quantization axis (for dipole mode)
    positions  :: Vector{Vec3}     # Positions of sites in global coordinates (Å)
    recipvecs  :: Mat3             # Reciprocal vectors in global coordinates (1/Å)
    energy_ϵ   :: Float64          # Energy shift applied to dynamical matrix prior to Bogoliubov transformation
    energy_tol :: Float64          # Energy tolerance for maximal imaginary part of quasiparticle energies

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    # dipole_corrs :: Bool                                  # Whether using all correlations from dipoles 
    # observables  :: Array{ComplexF64, 3}                  # Operators corresponding to observables
    # idxinfo      :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)
end

function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    # modename = swt.dipole_corrs ? "Dipole correlations" : "Custom correlations"
    modename = "Dipole correlations"
    printstyled(io, "SpinWaveTheory [$modename]\n"; bold=true, color=:underline)
    println(io, "Atoms in magnetic supercell: $(natoms(swt.sys.crystal))")
end

function num_bands(swt::SpinWaveTheory)
    (; sys) = swt
    nbosons = sys.mode == :SUN ? sys.Ns[1]-1 : 1
    return nbosons * natoms(sys.crystal)
end


"""
    dipole_to_angles

convert the dipole expectation values from the Cartesian frame to the spherical frame
"""
function dipole_to_angles(dipoles :: AbstractVector{Float64})
    r = norm(dipoles)
    @assert r > 1e-7
    θ = acos(dipoles[3] / r)
    @assert isfinite(θ)
    ϕ = atan(dipoles[2], dipoles[1])
    @assert isfinite(ϕ)
    (ϕ < 0.0) && (ϕ += 2π)
    return θ, ϕ
end

"""
    generate_local_sun_gens

Compute SU(N) generators in the local reference frame (for :SUN mode).
"""
# DD: Redo this using existing operator rotation facilities.
function generate_local_sun_gens(sys :: System)
    Nₘ, N = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert 
    S = (N-1)/2

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

    return s̃_mat, T̃_mat, Q̃_mat
end

"""
    generate_local_stevens_coefs

Compute the stevens coefficients in the local reference frame (for :dipole mode).
"""
function generate_local_stevens_coefs(sys :: System)
    c′_coef = Vector{StevensExpansion}()
    R_mat   = Vector{Mat3}()
    Nₘ = length(sys.dipoles) # number of magnetic atoms and dimension of Hilbert 
    R  = zeros(Float64, 3, 3)
    for atom = 1:Nₘ
        θ, ϕ = dipole_to_angles(sys.dipoles[1, 1, 1, atom])
        # U_mat = exp(-1im*ϕ*s_mat_N[3]) * exp(-1im*θ*s_mat_N[2])
        # SO(3) rotation that aligns the quantization axis. Important: since we will project out bosons that correspond to multipolar fluctuations,
        # therefore we use the explicit matrix to get rid of any ambiguity
        # Note that R * (0, 0, 1) = normalize(sys.dipoles[1,1,1,atom]))
        R[:] = [-sin(ϕ) -cos(ϕ)*cos(θ) cos(ϕ)*sin(θ);
                 cos(ϕ) -sin(ϕ)*cos(θ) sin(ϕ)*sin(θ);
                 0.0     sin(θ)        cos(θ)]

        (; c2, c4, c6) = sys.interactions_union[atom].onsite.stvexp

        SR  = Mat3(R)
        push!(R_mat, SR)
        c2′ = rotate_stevens_coefficients(c2, SR)
        c4′ = rotate_stevens_coefficients(c4, SR)
        c6′ = rotate_stevens_coefficients(c6, SR)
        c′  = StevensExpansion(c2′, c4′, c6′)
        push!(c′_coef, c′)
    end
    return R_mat, c′_coef
end


function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, energy_tol::Float64=1e-6) where N
    # Reshape into single unit cell
    cellsize_mag = cell_dimensions(sys) * diagm(collect(sys.latsize))
    sys = reshape_supercell_aux(sys, (1,1,1), cellsize_mag)

    # Computes the Stevens operator in the local reference frame and the SO(3) rotation matrix from global to local frame
    # (:dipole mode only)
    if sys.mode == :SUN
        s̃_mat, T̃_mat, Q̃_mat = generate_local_sun_gens(sys)
        c′_coef = Vector{StevensExpansion}()
        R_mat   = Vector{Mat3}()
    elseif sys.mode == :dipole
        R_mat, c′_coef = generate_local_stevens_coefs(sys)
        s̃_mat = zeros(ComplexF64, 0, 0, 0, 0)
        T̃_mat = zeros(ComplexF64, 0, 0, 0)
        Q̃_mat = zeros(ComplexF64, 0, 0, 0, 0)
    end

    positions = [global_position(sys, site) for site in all_sites(sys)][:]
    return SpinWaveTheory(sys, s̃_mat, T̃_mat, Q̃_mat, c′_coef, R_mat, positions, sys.crystal.recipvecs, energy_ϵ, energy_tol)
end
