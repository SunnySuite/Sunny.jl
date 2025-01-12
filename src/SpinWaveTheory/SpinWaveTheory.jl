struct SWTDataDipole
    local_rotations       :: Vector{Mat3}             # Rotations from global to quantization frame
    observables_localized :: Array{Vec3, 2}           # Observables rotated to local frame (nsites, nobs)
    stevens_coefs         :: Vector{StevensExpansion} # Rotated onsite coupling as Steven expansion
    sqrtS                 :: Vector{Float64}          # Square root of spin magnitudes
end

struct SWTDataSUN
    local_unitaries       :: Vector{Matrix{ComplexF64}} # Transformations from global to quantization frame
    observables_localized :: Array{HermitianC64, 2}     # Observables rotated to local frame (nobs × nsites)
    spins_localized       :: Array{HermitianC64, 2}     # Spins rotated to local frame (3 × nsites)
end

# To facilitate sharing some code with SpinWaveTheorySpiral
abstract type AbstractSpinWaveTheory end

"""
    SpinWaveTheory(sys::System; measure, regularization=1e-8)

Constructs an object to perform linear spin wave theory. The system must be in
an energy minimizing configuration. Enables calculation of [`dispersion`](@ref)
bands. If pair correlations are specified with `correspec`, one can also
calculate [`intensities_bands`](@ref) and broadened [`intensities`](@ref).

The spins in system must be energy-minimized, otherwise the Cholesky step of the
Bogoliubov diagonalization procedure will fail. The parameter `regularization`
adds a small positive shift to the diagonal of the dynamical matrix to avoid
numerical issues with quasi-particle modes of vanishing energy. Physically, this
shift can be interpreted as application of an inhomogeneous field aligned with
the magnetic ordering.
"""
struct SpinWaveTheory <: AbstractSpinWaveTheory
    sys            :: System
    data           :: Union{SWTDataDipole, SWTDataSUN}
    measure        :: MeasureSpec
    regularization :: Float64
end

function SpinWaveTheory(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8, energy_ϵ=nothing)
    if !isnothing(energy_ϵ)
        @warn "Keyword argument energy_ϵ is deprecated! Use `regularization` instead."
        regularization = energy_ϵ
    end

    measure = @something measure empty_measurespec(sys)
    if nsites(sys) != prod(size(measure.observables)[2:5])
        error("Size mismatch. Check that measure is built using consistent system.")
    end

    # Create single chemical cell that matches the full system size.
    new_shape = cell_shape(sys) * diagm(Vec3(sys.dims))
    new_cryst = reshape_crystal(orig_crystal(sys), new_shape)

    # Sort crystal positions so that their order matches sites in sys. Quadratic
    # scaling in system size.
    global_positions = global_position.(Ref(sys), vec(eachsite(sys)))
    p = map(new_cryst.positions) do r
        pos = new_cryst.latvecs * r
        findfirst(global_positions) do refpos
            isapprox(pos, refpos, atol=new_cryst.symprec)
        end
    end
    @assert allunique(p)
    permute_sites!(new_cryst, p)

    # Create a new system with dims (1,1,1). A clone happens in all cases.
    sys = reshape_supercell_aux(sys, new_cryst, (1,1,1))

    # Rotate local operators to quantization axis
    data = swt_data(sys, measure)

    return SpinWaveTheory(sys, data, measure, regularization)
end


function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    printstyled(io, "SpinWaveTheory ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", natoms(swt.sys.crystal), " atoms")
end

function nflavors(swt::SpinWaveTheory)
    (; sys) = swt
    nflavors = sys.mode == :SUN ? sys.Ns[1]-1 : 1
end

function nbands(swt::SpinWaveTheory)
    (; sys) = swt
    return nflavors(swt) * natoms(sys.crystal)
end


# Given q in reciprocal lattice units (RLU) for the original crystal, return a
# q_reshaped in RLU for the possibly-reshaped crystal.
function to_reshaped_rlu(sys::System{N}, q) where N
    return sys.crystal.recipvecs \ (orig_crystal(sys).recipvecs * q)
end

function dynamical_matrix(swt::SpinWaveTheory, q_reshaped)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    dynamical_matrix!(H, swt, q_reshaped)
    return Hermitian(H)
end

function dynamical_matrix!(H, swt::SpinWaveTheory, q_reshaped)
    if swt.sys.mode == :SUN
        swt_hamiltonian_SUN!(H, swt, q_reshaped)
    else
        @assert swt.sys.mode in (:dipole, :dipole_uncorrected)
        swt_hamiltonian_dipole!(H, swt, q_reshaped)
    end
end

function mul_dynamical_matrix(swt, x, qs_reshaped)
    y = zero(x)
    mul_dynamical_matrix!(swt, y, x, qs_reshaped)
end

function mul_dynamical_matrix!(swt, y, x, qs_reshaped)
    @assert size(x) == size(y) "Dimensions of x and y do not match"
    @assert size(x, 1) == length(qs_reshaped)

    if swt.sys.mode == :SUN
        multiply_by_hamiltonian_SUN!(y, x, swt, qs_reshaped)
    else
        multiply_by_hamiltonian_dipole!(y, x, swt, qs_reshaped)
    end

    return y
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
function swt_data(sys::System{N}, measure) where N
    # Calculate transformation matrices into local reference frames
    Na = nsites(sys)
    Nobs = size(measure.observables, 1)
    observables = reshape(measure.observables, Nobs, Na)

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = Vector{Matrix{ComplexF64}}(undef, Na)
    observables_localized = Array{HermitianC64}(undef, Nobs, Na)
    spins_localized = Array{HermitianC64}(undef, 3, Na)

    for i in 1:Na
        # Create unitary that rotates [0, ..., 0, 1] into ground state direction
        # Z that defines quantization axis
        Z = sys.coherents[i]

        U = if iszero(Z)
            # Set all operators on a vacant site to zero
            zeros(ComplexF64, N, N)
        else
            # Build unitary U satisfying U e_N ∝ Z.
            hcat(nullspace(Z'), Z)
        end
        local_unitaries[i] = U

        # Rotate observables into local reference frames
        for μ in 1:Nobs
            observables_localized[μ, i] = Hermitian(U' * observables[μ, i] * U)
        end

        S = spin_matrices_of_dim(; N)
        for μ in 1:3
            spins_localized[μ, i] = Hermitian(U' * S[μ] * U)
        end
    end

    # Rotate interactions into local reference frames
    for i in 1:Na
        Ui = local_unitaries[i]
        int = sys.interactions_union[i]

        # Zeeman coupling operator
        S = spin_matrices_of_dim(; N)
        B = sys.gs[i]' * sys.extfield[i]
        zeeman = B' * S

        # Merge and rotate all onsite couplings
        int.onsite = Hermitian(Ui' * (zeeman + int.onsite) * Ui)

        # Transform pair couplings into tensor decomposition and rotate.
        pair_new = PairCoupling[]
        for pc in int.pair
            # Convert PairCoupling to a purely general (tensor decomposed) interaction.
            pc_general = as_general_pair_coupling(pc, sys)

            # Rotate tensor decomposition into local frame.
            bond = pc.bond
            @assert bond.i == i
            Uj = local_unitaries[bond.j]
            pc_rotated = rotate_general_coupling_into_local_frame(pc_general, Ui, Uj)

            push!(pair_new, pc_rotated)
        end
        int.pair = pair_new
    end

    return SWTDataSUN(
        local_unitaries,
        observables_localized,
        spins_localized,
    )
end

# Compute Stevens coefficients in the local reference frame
function swt_data(sys::System{0}, measure)
    Na = nsites(sys)
    Nobs = size(measure.observables, 1)

    # Operators for rotating vectors into local frame
    Rs = map(1:Na) do i
        # Defines quantization axis.
        n = sys.dipoles[1, 1, 1, i]

        R = if iszero(n)
            # Set all operators on a vacant site to zero
            zero(Mat3)
        else
            # Build rotation R satisfying R z ∝ n.
            rotation_between_vectors([0, 0, 1], n)
        end

        # Rotation about the quantization axis is a U(1) gauge symmetry. The
        # angle θ below, for each atom, is arbitrary. We include this rotation
        # as an extra check of correctness.
        θ = 0.1 * i
        return R * axis_angle_to_matrix([0, 0, 1], θ)
    end

    # Operators for rotating Stevens quadrupoles into local frame
    Vs = map(Rs) do R
        iszero(R) ? zero(Mat5) : Mat5(operator_for_stevens_rotation(2, R))
    end

    # Observable is semantically a 1x3 row vector but stored in transpose
    # (column) form. To achieve effective right-multiplication by R, we should
    # in practice left-multiply column vector by R'.
    obs = reshape(measure.observables, Nobs, Na)
    obs_localized = [Rs[i]' * obs[μ, i] for μ in 1:Nobs, i in 1:Na]

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

    # Rotated Stevens expansion.
    cs = map(sys.interactions_union, Rs) do int, R
        iszero(R) ? empty_anisotropy(sys.mode, 0) : rotate_operator(int.onsite, R)
    end

    # Square root of spin magnitudes
    sqrtS = sqrt.(vec(sys.κs))

    return SWTDataDipole(Rs, obs_localized, cs, sqrtS)
end

# i is a site index for the flattened swt.sys. However, `measure` was originally
# constructed for a system with nontrivial lattice dims. Use fld1(i, prod(dims))
# to get a true atom index for the unflattened sys. This is needed to index
# measure.formfactors.
function get_swt_formfactor(measure, μ, i)
    sys_dims = size(measure.observables)[2:4]
    measure.formfactors[μ, fld1(i, prod(sys_dims))]
end
