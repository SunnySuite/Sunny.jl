###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################

"""
Additional fields for linear spin-wave calculations.
"""
struct SpinWaveTheory
    sys   :: System
    s̃_mat :: Array{ComplexF64, 4}  # dipole operators
    T̃_mat :: Array{ComplexF64, 3}  # single-ion anisos
    Q̃_mat :: Array{ComplexF64, 4}  # quarupolar operators (for biquad only)
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

Compute SU(N) generators in the local reference frame.
"""
# DD: Redo this using existing operator rotation facilities.
function generate_local_sun_gens(sys :: System)
    Nₘ, N = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    if sys.mode == :SUN
        s_mat = spin_matrices(N)

        s̃_mat = Array{ComplexF64, 4}(undef, N, N, 3, Nₘ)
        T̃_mat = Array{ComplexF64, 3}(undef, N, N, Nₘ)
        Q̃_mat = zeros(ComplexF64, 0, 0, 0, 0)

        U_mat = Matrix{ComplexF64}(undef, N, N)

        for atom = 1:Nₘ
            U_mat[:, 1] = sys.coherents[1, 1, 1, atom]
            U_mat[:, 2:N] = nullspace(U_mat[:, 1]')
            @assert isapprox(U_mat * U_mat', I) "rotation matrix from (global frame to local frame) not unitary"
            for μ = 1:3
                s̃_mat[:, :, μ, atom] = Hermitian(U_mat' * s_mat[μ] * U_mat)
            end
            T̃_mat[:, :, atom] = Hermitian(U_mat' * sys.interactions_union[atom].aniso.matrep * U_mat)
        end

    elseif sys.mode == :dipole
        s_mat_2 = spin_matrices(2)
        s_mat_N = spin_matrices(N)
        S = (N-1)/2

        # we support the biquad interactions now in the :dipole mode
        # we choose a particular basis of the nematic operators listed in Appendix B of *Phys. Rev. B 104, 104409*
        Q_mat = Vector{Matrix{ComplexF64}}(undef, 5)
        Q_mat[1] = -(s_mat_N[1] * s_mat_N[3] + s_mat_N[3] * s_mat_N[1])
        Q_mat[2] = -(s_mat_N[2] * s_mat_N[3] + s_mat_N[3] * s_mat_N[2])
        Q_mat[3] = s_mat_N[1] * s_mat_N[1] - s_mat_N[2] * s_mat_N[2]
        Q_mat[4] = s_mat_N[1] * s_mat_N[2] + s_mat_N[2] * s_mat_N[1]
        Q_mat[5] = √3 * s_mat_N[3] * s_mat_N[3] - 1/√3 * S * (S+1) * I
        
        s̃_mat = Array{ComplexF64, 4}(undef, 2, 2, 3, Nₘ)

        no_single_ion = isempty(sys.interactions_union[1].aniso.matrep)
        T̃_mat = no_single_ion ? zeros(ComplexF64, 0, 0, 0) : Array{ComplexF64, 3}(undef, 2, 2, Nₘ)
        Q̃_mat = Array{ComplexF64, 4}(undef, 2, 2, 5, Nₘ)

        U_mat_2 = Matrix{ComplexF64}(undef, 2, 2)
        U_mat_N = Matrix{ComplexF64}(undef, N, N)

        for atom = 1:Nₘ
            θ, ϕ = dipole_to_angles(sys.dipoles[1, 1, 1, atom])
            U_mat_N[:] = exp(-1im * ϕ * s_mat_N[3]) * exp(-1im * θ * s_mat_N[2])
            U_mat_2[:] = exp(-1im * ϕ * s_mat_2[3]) * exp(-1im * θ * s_mat_2[2])
            @assert isapprox(U_mat_N * U_mat_N', I) "rotation matrix from (global frame to local frame) not unitary"
            @assert isapprox(U_mat_2 * U_mat_2', I) "rotation matrix from (global frame to local frame) not unitary"
            for μ = 1:3
                s̃_mat[:, :, μ, atom] = Hermitian(U_mat_2' * s_mat_2[μ] * U_mat_2)
            end
            for ν = 1:5
                Q̃_mat[:, :, ν, atom] = Hermitian(U_mat_N' * Q_mat[ν] * U_mat_N)[1:2, 1:2]
            end

            if !no_single_ion
                T̃_mat[:, :, atom] = Hermitian(U_mat_N' * sys.interactions_union[atom].aniso.matrep * U_mat_N)[1:2, 1:2]
            end
        end
    end

    return s̃_mat, T̃_mat, Q̃_mat
end

"""
External constructor for `SpinWaveTheory`
"""
function SpinWaveTheory(sys::System, energy_ϵ::Float64=1e-8, energy_tol::Float64=1e-6)
    # Reshape into single unit cell
    cellsize_mag = cell_dimensions(sys) * diagm(collect(sys.latsize))
    sys = reshape_geometry_aux(sys, (1,1,1), cellsize_mag)

    s̃_mat, T̃_mat, Q̃_mat = generate_local_sun_gens(sys)

    latvecs_mag = isnothing(sys.origin) ? diagm(ones(3)) : sys.origin.crystal.latvecs \ sys.crystal.latvecs # DD: correct/necessary? 
    positions_chem = Vec3.([latvecs_mag * position for position in sys.crystal.positions]) # Positions of atoms in chemical coordinates
    recipvecs_mag = inv(latvecs_mag)'
    latvecs = isnothing(sys.origin) ? diagm(ones(3)) : sys.origin.crystal.latvecs # DD: correct/necessary?
    recipvecs_chem = inv(latvecs)'

    return SpinWaveTheory(sys, s̃_mat, T̃_mat, Q̃_mat, positions_chem, recipvecs_chem, recipvecs_mag, energy_ϵ, energy_tol)
end

"""
    chemical_to_magnetic

Convert the components of a wavevector from the original Brillouin zone (of the chemical lattice) to the reduced Brillouin zone (BZ)
(of the magnetic lattice). \
This is necessary because components in the reduced BZ are good quantum numbers.
`K` is the reciprocal lattice vector, and `k̃` is the components of wavevector in the reduced BZ. Note `k = K + k̃`
"""
function chemical_to_magnetic(sw_fields :: SpinWaveTheory, k)
    k = Vec3(k)
    α = sw_fields.recipvecs_mag \ k
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
    k_check = sw_fields.recipvecs_mag * (K + k̃)
    @assert norm(k - k_check) < 1e-12

    return K, k̃
end