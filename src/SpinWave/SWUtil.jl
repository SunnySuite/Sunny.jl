###########################################################################
# Below takes Sunny to construct `SpinWaveFields` for LSWT calculations.  #
###########################################################################

"""
Additional fields for linear spin-wave calculations.
"""
struct SpinWaveFields
    sys   :: Sunny.System
    s̃_mat :: Array{ComplexF64, 4}  # dipole operators
    T̃_mat :: Array{ComplexF64, 3}  # single-ion anisos
    Q̃_mat :: Array{ComplexF64, 4}  # quarupolar operators (for biquad only)
    chemical_positions :: Vector{Sunny.Vec3} # positions of magnetic atoms in units of (a₁, a₂, a₃) of the chemical lattice. (useful when computing the dynamical spin structure factor)
    chemic_reciprocal_basis :: Sunny.Mat3 # maybe not useful if we have David's interface for S(q, ω)
    maglat_reciprocal_basis :: Sunny.Mat3 # reciprocal lattice basis vectors for the magnetic supercell
    energy_ϵ   :: Float64 # energy epsilon in the diagonalization. Set to add to diagonal elements of the spin-wave Hamiltonian for cholesky decompostion
    energy_tol :: Float64 # energy tolerance for maximal imaginary part of spin-wave energies
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
function generate_local_sun_gens(sys :: Sunny.System)
    Nₘ, N = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    spin  = (N-1)/2
    if sys.mode == :SUN
        s_mat = Sunny.spin_matrices(N)

        s̃_mat = Array{ComplexF64, 4}(undef, N, N, 3, Nₘ)
        T̃_mat = Array{ComplexF64, 3}(undef, N, N, Nₘ)

        U_mat = Matrix{ComplexF64}(undef, N, N)

        for site = 1:Nₘ
            U_mat[:, 1] = sys.coherents[1, 1, 1, site]
            U_mat[:, 2:N] = nullspace(U_mat[:, 1]')
            @assert isapprox(U_mat * U_mat', I) "rotation matrix from (global frame to local frame) not unitary"
            for μ = 1:3
                s̃_mat[:, :, μ, site] = Hermitian(U_mat' * s_mat[μ] * U_mat)
            end
            T̃_mat[:, :, site] = Hermitian(U_mat' * sys.interactions.anisos[site].matrep * U_mat)
        end

    elseif sys.mode == :dipole
        s_mat_2 = spin * Sunny.spin_matrices(2)
        s_mat_N = Sunny.spin_matrices(N)
        
        s̃_mat = Array{ComplexF64, 4}(undef, 2, 2, 3, Nₘ)
        T̃_mat = Array{ComplexF64, 3}(undef, 2, 2, Nₘ)

        U_mat_2 = Matrix{ComplexF64}(undef, 2, 2)
        U_mat_N = Matrix{ComplexF64}(undef, N, N)

        for site = 1:Nₘ
            θ, ϕ = dipole_to_angles(sys.dipoles[1, 1, 1, site])
            U_mat_N[:] = exp(-1im * ϕ * s_mat_N[3]) * exp(-1im * θ * s_mat_N[2])
            U_mat_2[:] = exp(-1im * ϕ * s_mat_2[3]) * exp(-1im * θ * s_mat_2[2])
            @assert isapprox(U_mat_N * U_mat_N', I) "rotation matrix from (global frame to local frame) not unitary"
            @assert isapprox(U_mat_2 * U_mat_2', I) "rotation matrix from (global frame to local frame) not unitary"
            for μ = 1:3
                s̃_mat[:, :, μ, site] = Hermitian(U_mat_2' * s_mat_2[μ] * U_mat_2)
            end
            T̃_mat[:, :, site] = Hermitian(U_mat_N' * sys.interactions.anisos[site].matrep * U_mat_N)[1:2, 1:2]
        end
    
    # here I need help from Kipton and David to map back the solution, because in general this is difficult to implement.
    elseif sys.mode == :large_S
        error("unsupported")
        # spin  = (N-1) / 2
        # s_mat_2 = spin * Sunny.spin_matrices(2)

        # s̃_mat = Array{ComplexF64, 4}(undef, 2, 2, 3, Nₘ)
        # T̃_mat = Array{ComplexF64, 3}(undef, 2, 2, Nₘ)

        # U_mat_2 = Matrix{ComplexF64}(undef, 2, 2)

        # for site = 1:Nₘ
        #     θ, ϕ = dipole_to_angles(sys.dipoles[1, 1, 1, site])
        #     U_mat_2[:] = exp(-1im * ϕ * s_mat_2[3]) * exp(-1im * θ * s_mat_2[2])
        #     @assert isapprox(U_mat_2 * U_mat_2', I) "rotation matrix from (global frame to local frame) not unitary"
        #     for μ = 1:3
        #         s̃_mat[:, :, μ, site] = Hermitian(U_mat_2' * s_mat_2[μ] * U_mat_2)
        #     end
        #     # T̃_mat[:, :, site] = Hermitian(U_mat_N' * sys.interactions.anisos[site].matrep * U_mat_N)[1:2, 1:2]
        # end
    end

    return s̃_mat, T̃_mat
end

"""
External constructor for `SpinWaveFields`
"""
function SpinWaveFields(sys :: Sunny.System, energy_ϵ :: Float64=1e-8, energy_tol :: Float64=1e-6)
    s̃_mat, T̃_mat = generate_local_sun_gens(sys)
    Q̃_mat = zeros(ComplexF64, 0, 0, 0, 0)
    maglat_basis = sys.origin.crystal.lat_vecs \ sys.crystal.lat_vecs

    Nₘ = length(sys.dipoles)
    chemical_positions = Vector{Sunny.Vec3}(undef, Nₘ)
    for site = 1:Nₘ
        tmp_pos = maglat_basis * sys.crystal.positions[site]
        chemical_positions[site] = Sunny.Vec3(tmp_pos)
    end

    # computes the reciprocal basis vectors of the magnetic lattice (units 2π/|a|)
    det_A = det(maglat_basis')
    maglat_reciprocal_basis = zeros(Float64, 3, 3)
    maglat_reciprocal_basis[:, 1] = cross(maglat_basis[:, 2], maglat_basis[:, 3]) / det_A
    maglat_reciprocal_basis[:, 2] = cross(maglat_basis[:, 3], maglat_basis[:, 1]) / det_A
    maglat_reciprocal_basis[:, 3] = cross(maglat_basis[:, 1], maglat_basis[:, 2]) / det_A
    maglat_reciprocal_basis = Sunny.Mat3(maglat_reciprocal_basis)

    # computes the reciprocal basis vectors of the chemical lattice (units Å⁻¹)
    lat_vecs = sys.origin.crystal.lat_vecs
    det_A = det(lat_vecs')
    chemic_reciprocal_basis = zeros(Float64, 3, 3)
    chemic_reciprocal_basis[:, 1] = cross(lat_vecs[:, 2], lat_vecs[:, 3]) / det_A
    chemic_reciprocal_basis[:, 2] = cross(lat_vecs[:, 3], lat_vecs[:, 1]) / det_A
    chemic_reciprocal_basis[:, 3] = cross(lat_vecs[:, 1], lat_vecs[:, 2]) / det_A
    chemic_reciprocal_basis = Sunny.Mat3(chemic_reciprocal_basis)

    return SpinWaveFields(sys, s̃_mat, T̃_mat, Q̃_mat, chemical_positions, chemic_reciprocal_basis, maglat_reciprocal_basis, energy_ϵ, energy_tol)
end

"""
    k_chemical_to_k_magnetic

Convert the components of a wavevector from the original Brillouin zone (of the chemical lattice) to the reduced Brillouin zone (BZ)
(of the magnetic lattice). \
This is necessary because components in the reduced BZ are good quantum numbers.
`K` is the reciprocal lattice vector, and `k̃` is the components of wavevector in the reduced BZ. Note `k = K + k̃`
"""
function k_chemical_to_k_magnetic(sw_fields :: SpinWaveFields, k :: Vector{Float64})
    k_copy = deepcopy(k)
    α = sw_fields.maglat_reciprocal_basis \ k_copy
    k̃ = Vector{Float64}(undef, 3)
    K = Vector{Int}(undef, 3)
    for i = 1:3
        if abs(α[i]) < eps()
            K[i] = 0.0
            k̃[i] = 0.0
        else
            K[i] = Int(round(floor(α[i])))
            k̃[i] = α[i] - K[i]
        end
        @assert k̃[i] ≥ 0.0 && k̃[i] < 1.0
    end
    k_check = sw_fields.maglat_reciprocal_basis * (K + k̃)
    @assert norm(k - k_check) < 1e-12

    return K, k̃
end

function construct_magnetic_supercell(sys :: Sunny.System, A :: Matrix{Int})
    newsys = reshape_geometry(sys, A)
    mag_latsize = (1, 1, 1)
    mag_cell_size = Sunny.cell_dimensions(newsys) * diagm(collect(newsys.latsize))
    return Sunny.reshape_geometry_aux(newsys, mag_latsize, mag_cell_size)
end

"""
    set_dipoles!

Manually set the ground state configurations by inputing the `dipoles` for all sites along with their `positions` in the magnetic supercell. `dipoles` and `positions` are `3×Nₘ` matrices. Each column of the two objects should match each other. \n
*Warning*: 1. Must run `construct_magnetic_supercell` before calling this function. 2. You should always use the Langevin sampler provided by Sunny to find the ground state other than calling this function, unless the ground state configuration is known for sure. 3. This function works for `:dipole` and `large_S` mode. For `:SUN` mode, only the fundamental representation of `SU(2)`, i.e. S=1/2 works.
"""
function set_dipoles!(sys :: Sunny.System, dipoles :: Matrix{Float64}, positions :: Matrix{Float64})
    Ns, Nₘ = sys.Ns[1], length(sys.dipoles)
    spin = (Ns-1) / 2
    (sys.mode == :SUN && !isapprox(spin, 0.5, atol=1e-6)) && (throw("SU(N)N>2 not supported, use set_coherents! instead."))
    @assert sys.latsize == (1, 1, 1) "`sys` is not a valid magnetic supercell, run `construct_magnetic_supercell` first!"
    @assert size(dipoles, 1) == 3 "dipole should be a 3-vector"
    @assert size(positions, 1) == 3 "dipole should be a 3-vector"
    @assert size(dipoles, 2) == Nₘ "number of dipoles does not match `sys`"
    @assert size(positions, 2) == Nₘ "number of dipoles does not match `sys`"
    
    indices = Vector{Int}()

    for i in 1:Nₘ
        idx = findfirst(isapprox(positions[:, i], atol=1e-8, rtol=1e-10), sys.crystal.positions)
        @assert !isnothing(idx) "wrong positions for magnetic sites"
        push!(indices, idx)
    end

    for i in 1:Nₘ
        sys.dipoles[1, 1, 1, i] = dipoles[:, indices[i]]
    end

end

# """
#     set_coherents!

# Manually set the ground state configurations by inputing the `dipoles` for all sites along with their `positions` in the magnetic supercell. `dipoles` and `positions` are `3×Nₘ` matrices. Each column of the two objects should match each other. \n
# *Warning*: 1. Must run `construct_magnetic_supercell` before calling this function. 2. You should always use the Langevin sampler provided by Sunny to find the ground state other than calling this function, unless the ground state configuration is known for sure. 3. This function works for `:dipole` and `large_S` mode. For `:SUN` mode, only the fundamental representation of `SU(2)`, i.e. S=1/2 works.
# """
# function set_dipoles!(sys :: Sunny.System, dipoles :: Matrix{Float64}, positions :: Matrix{Float64})
#     Ns, Nₘ = sys.Ns[1], length(sys.dipoles)
#     spin = (Ns-1) / 2
#     (sys.mode == :SUN && !isapprox(spin, 0.5, atol=1e-6)) && (throw("SU(N)N>2 not supported, use set_coherents! instead."))
#     @assert sys.latsize == (1, 1, 1) "`sys` is not a valid magnetic supercell, run `construct_magnetic_supercell` first!"
#     @assert size(dipoles, 1) == 3 "dipole should be a 3-vector"
#     @assert size(positions, 1) == 3 "dipole should be a 3-vector"
#     @assert size(dipoles, 2) == Nₘ "number of dipoles does not match `sys`"
#     @assert size(positions, 2) == Nₘ "number of dipoles does not match `sys`"
    

#     indices = Vector{Int}(undef, Nₘ)

#     for i in 1:Nₘ
#         idx = findfirst(isapprox(positions[i], atol=1e-8, rtol=1e-10), sys.crystal.positions)
#         @assert !isnothing(idx) "wrong positions for magnetic sites"
#         push!(indices, idx)
#     end

#     sys.dipoles = dipoles[indices]

# end