struct SWTDataDipole
    local_rotations       :: Vector{Mat3}             # Rotations from global to local frame
    observables_localized :: Array{ComplexF64, 4}     # Observables in local frame
    stevens_coefs         :: Vector{StevensExpansion} # Rotated onsite coupling as Steven expansion
    sqrtS                 :: Vector{Float64}          # Square root of spin magnitudes
end

struct SWTDataSUN
    local_unitaries       :: Array{ComplexF64, 3}     # Aligns to quantization axis on each site
    spins_localized       :: Array{ComplexF64, 4}     # Spins in local frame
    observables_localized :: Array{ComplexF64, 4}     # Observables in local frame
end

"""
    SpinWaveTheory(sys, energy_ϵ::Float64=1e-8)

Constructs an object to perform linear spin wave theory. Use it with
[`dispersion`](@ref) and [`dssf`](@ref) functions.

The optional parameter `energy_ϵ` adds a small positive shift to the diagonal of
the dynamical matrix ``D`` to avoid numerical issues with zero-energy
quasi-particle modes.
"""
struct SpinWaveTheory
    sys          :: System
    data         :: Union{SWTDataDipole, SWTDataSUN}
    energy_ϵ     :: Float64
    observables  :: ObservableInfo
end

function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, observables=nothing, correlations=nothing, apply_g=true) where N
    # Reshape into single unit cell. A clone will always be performed, even if
    # no reshaping happens.
    cellsize_mag = cell_shape(sys) * diagm(Vec3(sys.latsize))
    sys = reshape_supercell_aux(sys, (1,1,1), cellsize_mag)

    # Rotate local operators to quantization axis
    if sys.mode == :SUN
        obs = parse_observables(N; observables, correlations, g = apply_g ? sys.gs : nothing)
        data = swt_data_sun(sys, obs)
    else
        if !isnothing(observables) || !isnothing(correlations)
            error("Only the default spin operators are supported in dipole mode")
        end
        obs = parse_observables(N; observables, correlations=nothing, g = apply_g ? sys.gs : nothing)
        data = swt_data_dipole(sys, obs)
    end

    return SpinWaveTheory(sys, data, energy_ϵ, obs)
end


function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    printstyled(io, "SpinWaveTheory\n"; bold=true, color=:underline)
    println(io, "Atoms in magnetic supercell: $(natoms(swt.sys.crystal))")
    show(io,MIME("text/plain"),swt.observables)
end

function nbands(swt::SpinWaveTheory)
    (; sys) = swt
    nflavors = sys.mode == :SUN ? sys.Ns[1]-1 : 1
    return nflavors * natoms(sys.crystal)
end


# Given q in reciprocal lattice units (RLU) for the original crystal, return a
# q_reshaped in RLU for the possibly-reshaped crystal.
function to_reshaped_rlu(sys::System{N}, q) where N
    return sys.crystal.recipvecs \ (orig_crystal(sys).recipvecs * q)
end

# Take PairCoupling `pc` and use it to make a new, equivalent PairCoupling that
# contains all information about the interaction in the `general` (tensor
# decomposition) field.
function as_general_pair_coupling(pc, sys)
    (; bond, scalar, bilin, biquad, general) = pc
    N1 = sys.Ns[bond.i]
    N2 = sys.Ns[bond.j]

    accum = zeros(ComplexF64, N1*N2, N1*N2)

    # Add scalar part
    accum += scalar * I

    # Add bilinear part
    S1, S2 = to_product_space(spin_matrices((N1-1)/2), spin_matrices((N2-1)/2))
    J = bilin isa Float64 ? bilin*I(3) : bilin
    accum += S1' * J * S2

    # Add biquadratic part
    K = biquad isa Float64 ? diagm(biquad * Sunny.scalar_biquad_metric) : biquad
    O1, O2 = to_product_space(stevens_matrices_of_dim(2; N=N1), stevens_matrices_of_dim(2; N=N2))
    accum += O1' * K * O2

    # Add general part
    for (A, B) in general.data
        accum += kron(A, B) 
    end

    # Generate new interaction with extract_parts=false 
    scalar, bilin, biquad, general = decompose_general_coupling(accum, N1, N2; extract_parts=false)

    return PairCoupling(bond, scalar, bilin, biquad, general)
end

function rotate_general_coupling_into_local_frame(pc, U1, U2)
    (; bond, scalar, bilin, biquad, general) = pc
    data_new = map(general.data) do (A, B)
        (Hermitian(U1'*A*U1), Hermitian(U2'*B*U2))
    end
    td = TensorDecomposition(general.gen1, general.gen2, data_new)
    return PairCoupling(bond, scalar, bilin, biquad, td)
end

# Prepare local operators and observables for SU(N) spin wave calculation by
# rotating these into the local reference frame determined by the ground state.
function swt_data_sun(sys::System{N}, obs) where N
    # Calculate transformation matrices into local reference frames
    n_magnetic_atoms = natoms(sys.crystal)

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = zeros(ComplexF64, N, N, n_magnetic_atoms)
    spins_localized = zeros(ComplexF64, N, N, 3, n_magnetic_atoms)
    observables_localized = zeros(ComplexF64, N, N, num_observables(obs), n_magnetic_atoms)

    for atom in 1:n_magnetic_atoms
        # Create unitary that rotates [0, ..., 0, 1] into ground state direction
        # Z that defines quantization axis
        Z = sys.coherents[atom]
        view(local_unitaries, :, N, atom)     .= Z
        view(local_unitaries, :, 1:N-1, atom) .= nullspace(Z')
    end

    for atom in 1:n_magnetic_atoms
        U = view(local_unitaries, :, :, atom)

        # Rotate observables into local reference frames
        S = spin_matrices_of_dim(; N)
        for k in 1:3
            spins_localized[:, :, k, atom] = Hermitian(U' * S[k] * U)
        end
        
        for k in 1:num_observables(obs)
            A = observable_at_site(obs.observables[k],CartesianIndex(1,1,1,atom))
            observables_localized[:, :, k, atom] = Hermitian(U' * convert(Matrix, A) * U)
        end

        # Rotate interactions into local reference frames
        int = sys.interactions_union[atom]

        # Accumulate Zeeman terms into OnsiteCoupling
        S = spin_matrices_of_dim(; N)
        B = sys.units.μB * (sys.gs[atom]' * sys.extfield[atom])
        int.onsite += B' * S

        # Rotate onsite anisotropy
        int.onsite = Hermitian(U' * int.onsite * U) 

        # Transform pair couplings into tensor decomposition and rotate.
        pair_new = PairCoupling[]
        for pc in int.pair
            # Convert PairCoupling to a purely general (tensor decomposed) interaction.
            pc_general = as_general_pair_coupling(pc, sys)

            # Rotate tensor decomposition into local frame.
            bond = pc.bond
            @assert bond.i == atom
            U′ = view(local_unitaries, :, :, bond.j)
            pc_rotated = rotate_general_coupling_into_local_frame(pc_general, U, U′)

            push!(pair_new, pc_rotated)
        end
        int.pair = pair_new
    end

    return SWTDataSUN(
        local_unitaries,
        spins_localized,
        observables_localized
    )
end


# Compute Stevens coefficients in the local reference frame
function swt_data_dipole(sys::System{0}, obs)
    cs = StevensExpansion[]
    Rs = Mat3[]
    Vs = Mat5[]

    n_magnetic_atoms = natoms(sys.crystal)
    observables_localized = zeros(ComplexF64, 1, 3, num_observables(obs), n_magnetic_atoms)

    for atom in 1:n_magnetic_atoms
        # Direction n of dipole will define rotation R that aligns the
        # quantization axis.
        n = normalize(sys.dipoles[1,1,1,atom])

        # Build matrix that rotates from z to n.
        R = rotation_between_vectors([0, 0, 1], n)
        @assert R * [0, 0, 1] ≈ n

        # Rotation about the quantization axis is a U(1) gauge symmetry. The
        # angle θ below, for each atom, is arbitrary. We include this rotation
        # as an extra check of correctness.
        θ = 0.1 * atom
        R = R * axis_angle_to_matrix([0, 0, 1], θ)

        # Rotated Stevens expansion.
        c = rotate_operator(sys.interactions_union[atom].onsite, R)

        # Precalculate operators for rotating Stevens coefficients.
        V = operator_for_stevens_rotation(2, R)

        # Observables are 1x3 row vectors
        for μ = 1:num_observables(obs)
            row = observable_at_site(obs.observables[μ],CartesianIndex(1,1,1,atom))
            observables_localized[:,:,μ,atom] .= Matrix(row) * R
        end

        push!(Rs, R)
        push!(cs, c)
        push!(Vs, V)
    end

    # Precompute transformed exchange matrices and store in sys.interactions_union.
    for ints in sys.interactions_union
        for c in eachindex(ints.pair)
            (; bond, scalar, bilin, biquad, general) = ints.pair[c]
            (; i, j) = bond

            if !iszero(bilin)  # Leave zero if already zero
                J = Mat3(bilin*I)
                bilin = Rs[i]' * J * Rs[j]
            end

            if !iszero(biquad)
                J = biquad
                J = Mat5(J isa Number ? J * diagm(scalar_biquad_metric) : J)
                biquad = Vs[i]' * J * Vs[j]
            end

            @assert isempty(general.data)

            ints.pair[c] = PairCoupling(bond, scalar, bilin, biquad, general)
        end
    end

    sqrtS = [sqrt((N-1)/2) for N in vec(sys.Ns)]
    return SWTDataDipole(Rs, observables_localized, cs, sqrtS)
end
