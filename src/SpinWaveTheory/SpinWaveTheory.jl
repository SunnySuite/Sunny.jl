struct SWTDataDipole
    local_rotations       :: Vector{Mat3}
    stevens_coefs         :: Vector{StevensExpansion}
end

struct SWTDataSUN
    local_unitaries       :: Array{ComplexF64, 3}  # Aligns to quantization axis on each site
    zeeman_operators      :: Array{ComplexF64, 3}  # Consider constructing on the fly, as for dipole mode
    observable_operators  :: Array{ComplexF64, 4}  # Observables in local frame (for intensity calcs)
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

function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, observables=nothing, correlations=nothing) where N
    if !isnothing(sys.ewald)
        error("SpinWaveTheory does not yet support long-range dipole-dipole interactions.")
    end

    # Reshape into single unit cell. A clone will always be performed, even if
    # no reshaping happens.
    cellsize_mag = cell_shape(sys) * diagm(Vec3(sys.latsize))
    sys = reshape_supercell_aux(sys, (1,1,1), cellsize_mag)

    # Rotate local operators to quantization axis
    if sys.mode == :SUN
        obs = parse_observables(N; observables, correlations)
        data = swt_data_sun(sys, obs)
    else
        if !isnothing(observables) || !isnothing(correlations)
            error("Only the default spin operators are supported in dipole mode")
        end
        obs = parse_observables(N; observables, correlations=nothing)
        data = swt_data_dipole!(sys)
    end

    return SpinWaveTheory(sys, data, energy_ϵ, obs)
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


# Convert 3-vector from the Cartesian frame to the spherical frame.
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

# Take PairCoupling `pc` and use it to make a new, equivalent PairCoupling that
# contains all information about the interaction in the `general` (tensor
# decomposition) field.
function as_general_pair_coupling(pc, sys)
    (; isculled, bond, scalar, bilin, biquad, general) = pc
    N1 = sys.Ns[1, 1, 1, bond.i]
    N2 = sys.Ns[1, 1, 1, bond.j]

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

    return PairCoupling(isculled, bond, scalar, bilin, biquad, general)
end

function rotate_general_coupling_into_local_frame(pc, U1, U2)
    (; isculled, bond, scalar, bilin, biquad, general) = pc
    data_new = Tuple{HermitianC64, HermitianC64}[]
    for (A, B) in general.data
        push!(data_new, (Hermitian(U1'*A*U1), Hermitian(U2'*B*U2)))
    end
    td = TensorDecomposition(general.gen1, general.gen2, data_new)
    return PairCoupling(isculled, bond, scalar, bilin, biquad, td)
end

# Prepare local operators and observables for SU(N) spin wave calculation by
# rotating these into the local reference frame determined by the ground state.
function swt_data_sun(sys::System{N}, obs) where N
    # Calculate transformation matrices into local reference frames
    n_magnetic_atoms = natoms(sys.crystal)

    local_unitaries = Array{ComplexF64}(undef, N, N, n_magnetic_atoms)
    for atom in 1:n_magnetic_atoms
        basis = view(local_unitaries, :, :, atom)
        # First axis of local quantization basis is along the 
        # ground-state polarization axis
        basis[:, N] .= sys.coherents[1, 1, 1, atom]
        
        # Remaining axes are arbitrary but mutually orthogonal
        # and orthogonal to the first axis
        basis[:, 1:N-1] .= nullspace(basis[:, N]')
    end

    # Preallocate buffers for rotate operators and observables.
    zeeman_operators_localized = zeros(ComplexF64, N, N, n_magnetic_atoms)
    observables_localized = zeros(ComplexF64, N, N, num_observables(obs), n_magnetic_atoms)

    # Rotate observables into local reference frames and store (for intensities calculations).
    for atom in 1:n_magnetic_atoms
        U = view(local_unitaries, :, :,atom)
        for k = 1:num_observables(obs)
            observables_localized[:, :, k, atom] = Hermitian(U' * convert(Matrix,obs.observables[k]) * U)
        end
    end

    # Calculate external magnetic field operator on each site using the local
    # frames and accumulate into onsite_operator_localized.
    (; extfield, gs, units) = sys
    for atom in 1:n_magnetic_atoms
        N1 = sys.Ns[1, 1, 1, atom]
        B = units.μB * (gs[1, 1, 1, atom]' * extfield[1, 1, 1, atom])
        U = view(local_unitaries, :, :, atom) 
        Sx, Sy, Sz = map(Sa -> U' * Sa * U, spin_matrices_of_dim(; N=N1))
        @. zeeman_operators_localized[:, :, atom] -= B[1]*Sx + B[2]*Sy + B[3]*Sz
    end

    # Transform interactions inside the system into local reference frames,
    # and transform pair interactions into tensor decompositions.
    for idx in CartesianIndices(sys.interactions_union)
        atom = idx.I[end]   # Get index for unit cell, regardless of type of interactions_union
        int = sys.interactions_union[idx]

        # Rotate onsite anisotropy. (NB: could add Zeeman term into onsite, but
        # that might be confusing for maintenance/debugging and the performance
        # gains are minimal.)
        U = local_unitaries[:, :, atom]
        int.onsite = Hermitian(U' * int.onsite * U) 

        # Transform pair couplings into tensor decomposition and rotate.
        pair_new = PairCoupling[]
        for pc in int.pair
            # Convert PairCoupling to a purely general (tensor decomposed) interaction.
            pc_general = as_general_pair_coupling(pc, sys)

            # Rotate tensor decomposition into local frame.
            bond = pc.bond
            Ui, Uj = view(local_unitaries, :, :, bond.i), view(local_unitaries, :, :, bond.j)
            pc_rotated = rotate_general_coupling_into_local_frame(pc_general, Ui, Uj)

            push!(pair_new, pc_rotated)
        end
        int.pair = pair_new
    end

    return SWTDataSUN(
        local_unitaries,
        zeeman_operators_localized,
        observables_localized
    )
end


# Compute Stevens coefficients in the local reference frame
function swt_data_dipole!(sys::System{0})
    N = sys.Ns[1]
    S = (N-1)/2

    cs = StevensExpansion[]
    Rs = Mat3[]
    Vs = Mat5[]

    for atom in 1:natoms(sys.crystal)
        # SO(3) rotation that aligns the quantization axis. Important: since we
        # will project out bosons that correspond to multipolar fluctuations,
        # therefore we use the explicit matrix to get rid of any ambiguity.
        
        # As a unitary, U = exp(-i ϕ Sz) exp(-i θ Sy).
        θ, ϕ = dipole_to_angles(sys.dipoles[1,1,1,atom])
        R = SA[-sin(ϕ) -cos(ϕ)*cos(θ) cos(ϕ)*sin(θ);
                cos(ϕ) -sin(ϕ)*cos(θ) sin(ϕ)*sin(θ);
                0.0     sin(θ)        cos(θ)]

        # Rotated Stevens expansion.
        c = rotate_operator(sys.interactions_union[atom].onsite, R)

        # Precalculate operators for rotating Stevens coefficients.
        V = operator_for_stevens_rotation(2, R)

        push!(Rs, R)
        push!(cs, c)
        push!(Vs, V)
    end

    # Precompute transformed exchange matrices and store in sys.interactions_union.
    for ints in sys.interactions_union
        for c in eachindex(ints.pair)
            (; isculled, bond, scalar, bilin, biquad, general) = ints.pair[c]
            isculled && break
            i, j = bond.i, bond.j

            if !iszero(bilin)  # Leave zero if already zero
                J = Mat3(bilin*I)
                bilin = S * (Rs[i]' * J * Rs[j]) 
            end

            if !iszero(biquad)
                J = biquad
                J = Mat5(J isa Number ? J * diagm(scalar_biquad_metric) : J)
                biquad = S^3 * Mat5(Vs[i]' * J * Vs[j]) 
            end

            @assert isempty(general.data)

            ints.pair[c] = PairCoupling(isculled, bond, scalar, bilin, biquad, general)
        end
    end

    return SWTDataDipole(Rs, cs)
end
