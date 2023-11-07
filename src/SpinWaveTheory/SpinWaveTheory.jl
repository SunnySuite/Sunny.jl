###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################
struct SWTDataDipole
    local_rotations  :: Vector{Mat3}
    stevens_coefs    :: Vector{StevensExpansion}
end

struct SWTDataSUN
    dipole_operators      :: Array{ComplexF64, 4}
    quadrupole_operators  :: Array{ComplexF64, 4}
    onsite_operator       :: Array{ComplexF64, 3}  
    bond_operator_pairs   :: Vector{Tuple{Bond, Array{ComplexF64, 4}}}
    observable_operators  :: Array{ComplexF64, 4}
    local_unitary         :: Array{ComplexF64,3} # Aligns quantization axis on each site
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
    sys        :: System
    data       :: Union{SWTDataDipole, SWTDataSUN}
    energy_ϵ   :: Float64

    observables :: ObservableInfo
end

function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, observables=nothing, correlations=nothing) where N
    if !isnothing(sys.ewald)
        error("SpinWaveTheory does not yet support long-range dipole-dipole interactions.")
    end

    # Reshape into single unit cell. A clone will always be performed, even if
    # no reshaping happens.
    cellsize_mag = cell_shape(sys) * diagm(SVector(sys.latsize))
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

# Prepare local operators and observables for SU(N) spin wave calculation by
# rotating these into the local reference frame determined by the ground state.
function swt_data_sun(sys::System{N}, obs) where N
    # Calculate transformation matrices into local reference frames
    n_magnetic_atoms = natoms(sys.crystal)

    local_unitary = Array{ComplexF64}(undef, N, N, n_magnetic_atoms)
    for atom in 1:n_magnetic_atoms
        basis = view(local_unitary, :, :, atom)
        # First axis of local quantization basis is along the 
        # ground-state polarization axis
        basis[:, 1] .= sys.coherents[1, 1, 1, atom]
        
        # Remaining axes are arbitrary but mutually orthogonal
        # and orthogonal to the first axis
        basis[:, 2:N] .= nullspace(basis[:, 1]')
    end

    # Preallocate buffers for rotate operators and observables.
    dipole_operators_localized = zeros(ComplexF64, 3, N, N, n_magnetic_atoms)
    onsite_operator_localized = zeros(ComplexF64, N, N, n_magnetic_atoms)
    quadrupole_operators_localized = zeros(ComplexF64, 5, N, N, n_magnetic_atoms)
    observables_localized = zeros(ComplexF64, N, N, num_observables(obs), n_magnetic_atoms)

    # Rotate SU(N) bases and observables and store in dense array. Note that the
    # first index is the component index. As a result, (Sˣᵢ[m,n], Sʸᵢ[m,n], Sᶻᵢ[m,n]) is
    # stored contiguously for each matrix element mn. This is the natural order
    # for the current way of constructing the spin wave Hamiltonian.
    dipole_operators = spin_matrices_of_dim(; N)
    quadrupole_operators = stevens_matrices_of_dim(2; N)
    for atom in 1:n_magnetic_atoms
        U = view(local_unitary, :, :,atom)
        for μ = 1:3
            dipole_operators_localized[μ, :, :, atom] = Hermitian(U' * dipole_operators[μ] * U)
        end
        for ν = 1:5
            quadrupole_operators_localized[ν, :, :, atom] = Hermitian(U' * quadrupole_operators[ν] * U)
        end
        onsite_operator_localized[:, :, atom] = Hermitian(U' * sys.interactions_union[atom].onsite * U)
        for k = 1:num_observables(obs)
            observables_localized[:, :, k, atom] = Hermitian(U' * convert(Matrix,obs.observables[k]) * U)
        end
    end

    # Calculate external magnetic field operator on each site using the local
    # frames and accumulate into onsite_operator_localized.
    (; extfield, gs, units) = sys
    for atom in 1:n_magnetic_atoms
        B = units.μB * (gs[1, 1, 1, atom]' * extfield[1, 1, 1, atom])
        S = view(dipole_operators_localized, :, :, :, atom)
        @. onsite_operator_localized[:, :, atom] -= B[1]*S[1, :, :] + B[2]*S[2, :, :] + B[3]*S[3, :, :]
    end

    # Rotate operators generated by the tensor decomposition of generalized
    # interactions.
    bond_operator_pairs_localized = Tuple{Bond, Array{ComplexF64, 4}}[]
    for int in sys.interactions_union
        for pc in int.pair
            (; isculled, bond, general) = pc
            (isculled || length(general.data) == 0) && continue 

            Ui, Uj = local_unitary[:,:,bond.i], local_unitary[:,:,bond.j] 
            nops = length(general.data)
            bond_operators = zeros(ComplexF64, nops, N, N, 2)
            for (n, (A, B)) in enumerate(general.data)
                bond_operators[n,:,:,1] = Ui' * A * Ui 
                bond_operators[n,:,:,2] = Uj' * B * Uj 
            end
            push!(bond_operator_pairs_localized, (bond, bond_operators))
        end
    end

    return SWTDataSUN(
        dipole_operators_localized, 
        quadrupole_operators_localized, 
        onsite_operator_localized,
        bond_operator_pairs_localized,
        observables_localized,
        local_unitary
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
    for (n, ints) in enumerate(sys.interactions_union)
        for (c, coupling) in enumerate(ints.pair)
            (; isculled, bond, scalar, bilin, biquad, general) = coupling
            isculled && break
            i, j = bond.i, bond.j

            if !iszero(coupling.bilin)  # Leave zero if already zero
                J = Mat3(coupling.bilin*I)
                bilin = S * (Rs[i]' * J * Rs[j]) 
                sys.interactions_union[n].pair[c] = PairCoupling(isculled, bond, scalar, bilin, biquad, general)
            end

            if !iszero(coupling.biquad)
                J = coupling.biquad
                J = Mat5(J isa Number ? J * diagm(scalar_biquad_metric) : J)
                biquad = S^3 * Mat5(Vs[i]' * J * Vs[j]) 
                sys.interactions_union[n].pair[c] = PairCoupling(isculled, bond, scalar, bilin, biquad, general)
            end
        end
    end

    return SWTDataDipole(Rs, cs)
end