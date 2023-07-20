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
    s̃_mat :: Array{ComplexF64, 4}  # dipole operators
    T̃_mat :: Array{ComplexF64, 3}  # single-ion anisos
    Q̃_mat :: Array{ComplexF64, 4}  # quarupolar operators
    c′_coef :: Vector{StevensExpansion}    # c′_coefficents (for dipole mode)
    R_mat   :: Vector{Mat3}        # SO(3) rotation to align the quantization axis (for dipole mode)
    positions_chem :: Vector{Vec3} # positions of magnetic atoms in units of (a₁, a₂, a₃) of the chemical lattice. (useful when computing the dynamical spin structure factor)
    recipvecs_chem :: Mat3 # maybe not useful if we have David's interface for S(q, ω)
    recipvecs_mag :: Mat3 # reciprocal lattice basis vectors for the magnetic supercell
    energy_ϵ   :: Float64 # energy epsilon in the diagonalization. Set to add to diagonal elements of the spin-wave Hamiltonian for cholesky decompostion
    energy_tol :: Float64 # energy tolerance for maximal imaginary part of spin-wave energies

    # Correlation info (αβ indices of 𝒮^{αβ}(q,ω))
    # dipole_corrs :: Bool                                  # Whether using all correlations from dipoles 
    # observables  :: Array{ComplexF64, 3}                  # Operators corresponding to observables
    # idxinfo      :: SortedDict{CartesianIndex{2}, Int64}  # (α, β) to save from 𝒮^{αβ}(q, ω)
end

function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    # modename = swt.dipole_corrs ? "Dipole correlations" : "Custom correlations"
    modename = "Dipole correlations"
    printstyled(io, "SpinWaveTheory [$modename]\n"; bold=true, color=:underline)
    println(io, "Atoms in magnetic supercell $(length(swt.positions_chem))")
end

function num_bands(swt::SpinWaveTheory)
    (; sys) = swt
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nbands = Nf * Nm
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
    (ϕ < 0.0) && (ϕ += 2.0 * π)
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
        T̃_mat[:, :, atom] = Hermitian(U_mat' * sys.interactions_union[atom].aniso.matrep * U_mat)
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
    Nₘ, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert 
    s_mat_N = spin_matrices(N=Ns)
    Nₘ = length(sys.dipoles)
    s̃_mat = Array{ComplexF64, 4}(undef, 2, 2, 3, Nₘ)
    Nₘ = length(sys.dipoles)
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
        # In Cristian's note, S̃ = R S, so here we should pass SR'
        push!(R_mat, SR')
        for μ = 1:3
            s̃_mat[:, :, μ, atom] = Hermitian(rotate_operator(s_mat_N[μ], SR))[1:2, 1:2]
        end
        c2′ = rotate_stevens_coefficients(c2, SR)
        c4′ = rotate_stevens_coefficients(c4, SR)
        c6′ = rotate_stevens_coefficients(c6, SR)
        c′  = StevensExpansion(c2′, c4′, c6′)
        push!(c′_coef, c′)
    end
    return s̃_mat, R_mat, c′_coef
end


function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, energy_tol::Float64=1e-6) where N
    # Reshape into single unit cell
    cellsize_mag = cell_dimensions(sys) * diagm(collect(sys.latsize))
    sys = reshape_geometry_aux(sys, (1,1,1), cellsize_mag)

    # Computes the Stevens operator in the local reference frame and the SO(3) rotation matrix from global to local frame
    # (:dipole mode only)
    if sys.mode == :SUN
        s̃_mat, T̃_mat, Q̃_mat = generate_local_sun_gens(sys)
        c′_coef = Vector{StevensExpansion}()
        R_mat   = Vector{Mat3}()
    elseif sys.mode == :dipole
        s̃_mat, R_mat, c′_coef = generate_local_stevens_coefs(sys)
        T̃_mat = zeros(ComplexF64, 0, 0, 0)
        Q̃_mat = zeros(ComplexF64, 0, 0, 0, 0)
    end

    latvecs_mag = isnothing(sys.origin) ? diagm(ones(3)) : sys.origin.crystal.latvecs \ sys.crystal.latvecs # DD: correct/necessary? 
    positions_chem = Vec3.([latvecs_mag * position for position in sys.crystal.positions]) # Positions of atoms in chemical coordinates
    recipvecs_mag = inv(latvecs_mag)'
    latvecs_chem = isnothing(sys.origin) ? diagm(ones(3)) : sys.origin.crystal.latvecs # DD: correct/necessary?
    recipvecs_chem = inv(latvecs_chem)'

    return SpinWaveTheory(sys, s̃_mat, T̃_mat, Q̃_mat, c′_coef, R_mat, positions_chem, recipvecs_chem, recipvecs_mag, energy_ϵ, energy_tol)
end

"""
    chemical_to_magnetic

Convert the components of a wavevector from the original Brillouin zone (of the chemical lattice) to the reduced Brillouin zone (BZ)
(of the magnetic lattice). \
This is necessary because components in the reduced BZ are good quantum numbers.
`K` is the reciprocal lattice vector, and `k̃` is the components of wavevector in the reduced BZ. Note `k = K + k̃`
"""
function chemical_to_magnetic(swt::SpinWaveTheory, k)
    k = Vec3(k)
    α = swt.recipvecs_mag \ k
    k̃ = Vector{Float64}(undef, 3)
    K = Vector{Int}(undef, 3)
    for i = 1:3
        if abs(α[i]) < eps()
            K[i] = k̃[i] = 0.0
        else
            K[i] = Int(round(floor(α[i])))
            k̃[i] = α[i] - K[i]
        end
        @assert k̃[i] ≥ 0.0 && k̃[i] < 1.0
    end
    k_check = swt.recipvecs_mag * (K + k̃)
    @assert norm(k - k_check) < 1e-12

    return K, k̃
end
