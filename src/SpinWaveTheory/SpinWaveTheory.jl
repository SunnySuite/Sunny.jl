struct SWTDataDipole
    local_rotations       :: Vector{Mat3}             # Rotations from global to quantization frame
    observables_localized :: Array{Vec3, 2}           # Observables rotated to local frame (nobs × nsites)
    stevens_coefs         :: Vector{StevensExpansion} # Rotated onsite coupling as Steven expansion
    sqrtS                 :: Vector{Float64}          # Square root of spin magnitudes
end

struct SWTDataSUN
    local_unitaries  :: Vector{Matrix{ComplexF64}}   # Transformations from global to quantization frame
    obs_parts        :: Array{HermitianC64, 3}       # Rotated observable parts (ncontribs × nobs × nsites)
    observable_buf   :: Matrix{ComplexF64}           # Scratch buffer (N × N)
    spins_localized  :: Array{HermitianC64, 3}       # Spins in local frame (3 × nparts × nsites)
end

# To facilitate sharing some code with SpinWaveTheorySpiral
abstract type AbstractSpinWaveTheory end

"""
    SpinWaveTheory(sys::System; measure, regularization=1e-8)

Constructs an object to perform linear spin wave theory. The `measure` object
specifies observable fields, ``𝐪``-dependent contractions, and form factors.
Common choices are [`ssf_perp`](@ref), [`ssf_trace`](@ref), and
[`ssf_custom`](@ref). The resulting `SpinWaveTheory` object can be used to
calculate [`intensities_bands`](@ref) and broadened [`intensities`](@ref). If
`measure=nothing`, it is still possible to calculate the [`dispersion`](@ref)
curves. Eigenvectors describing the Bogoliubov bosons are available in
[`excitations`](@ref).

The magnetic structure in `sys` must be energy minimized, otherwise the Cholesky
step of the Bogoliubov diagonalization procedure will fail. The parameter
`regularization` adds a small positive shift to the diagonal of the dynamical
matrix to avoid numerical issues with quasi-particle modes of vanishing energy.
Physically, this shift can be interpreted as application of an inhomogeneous
field aligned with the magnetic ordering.
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
    if size(eachsite(sys)) != size(measure.operators)[3:6]
        error("Size mismatch. Check that measure is built using consistent system.")
    end

    # Create a new system with dims (1,1,1). A clone happens in all cases. For an
    # entangled system, `reshape_supercell_aux` flattens the physical
    # (uncontracted) system in tandem and rebuilds the entanglement metadata.
    new_cryst = resize_and_flatten_crystal(sys.crystal, sys.dims)
    sys = reshape_supercell_aux(sys, new_cryst, (1, 1, 1))

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


function to_standard_rlu(sys::System, q_reshaped)
    return orig_crystal(sys).recipvecs \ (sys.crystal.recipvecs * q_reshaped)
end

# Given q in reciprocal lattice units (RLU) for the original crystal, return a
# q_reshaped in RLU for the possibly-reshaped crystal.
function to_reshaped_rlu(sys::System, q)
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
    nparts = size(measure.operators, 1)
    Nobs = num_observables(measure)
    flat_ops = reshape(measure.operators, nparts, Nobs, Na)

    # Global-frame spin dipole "parts" for each site: the physical magnetic atoms
    # contributing to it, as `(spin_ops::NTuple{3, HermitianC64}, g::Mat3,
    # B::Vec3)`. For an ordinary system each site has one part (its bare spin S);
    # for an entangled unit the parts are the unit's bare atoms, each embedded into
    # the product space, with per-atom g-factor and field from the bare system.
    site_parts = swt_spin_parts(sys)
    nparts_spin = maximum(length, site_parts)

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = Vector{Matrix{ComplexF64}}(undef, Na)
    obs_parts = Array{HermitianC64}(undef, nparts, Nobs, Na)
    observable_buf = zeros(ComplexF64, N, N)
    spins_localized = fill(Hermitian(zeros(ComplexF64, N, N)), 3, nparts_spin, Na)

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
        for k in 1:nparts, μ in 1:Nobs
            obs_parts[k, μ, i] = Hermitian(U' * flat_ops[k, μ, i] * U)
        end

        # Spin dipole operators per part in the local frame (for dipole-dipole).
        for (k, (spin_ops, _, _)) in enumerate(site_parts[i])
            for α in 1:3
                spins_localized[α, k, i] = Hermitian(U' * spin_ops[α] * U)
            end
        end
    end

    # Rotate interactions into local reference frames
    for i in 1:Na
        Ui = local_unitaries[i]
        int = sys.interactions_union[i]

        # Zeeman coupling Σ_parts Σ_α (gₐ'Bₐ)_α Sₐᵅ, accumulated into the onsite
        # coupling in the global frame so it rotates with the anisotropy below.
        onsite = Matrix{ComplexF64}(int.onsite)
        for (spin_ops, g, B) in site_parts[i]
            gB = g' * B
            for α in 1:3
                onsite += gB[α] * spin_ops[α]
            end
        end
        int.onsite = Hermitian(Ui' * onsite * Ui)

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
        obs_parts,
        observable_buf,
        spins_localized,
    )
end

# The physical magnetic-atom "parts" contributing to each SU(N) site, as a list
# per site of `(spin_ops, g, B)` in the global frame. For an entangled system,
# the iterator runs over all sites of the bare system.
function swt_spin_parts(sys::System{N}) where N
    if is_entangled(sys)
        (; uncontracted, unit_map, bare_dipole_operators) = get_entanglement(sys)
        return map(unit_map.unit_to_members) do members
            map(members) do member
                (bare_dipole_operators[member.atom],
                 uncontracted.gs[1, 1, 1, member.atom],
                 uncontracted.extfield[1, 1, 1, member.atom])
            end
        end
    else
        return map(1:nsites(sys)) do i
            S = spin_matrices_of_dim(; N)
            [(ntuple(α -> S[α], 3), sys.gs[i], sys.extfield[i])]
        end
    end
end

# The system carrying physical magnetic moments for the long-range dipole-dipole
# term, and a map from each of its atoms to a `(site, part)` pair indexing
# `spins_localized`. For an ordinary system this is `sys` itself with the
# trivial map `a -> (a, 1)`. For an entangled system it is the physical bare
# system, with `unit_map.atom_to_unit[a] = (unit, part)` mapping each bare atom
# to its unit and intra-unit part index. Both share `sys.crystal.latvecs`, so
# the Ewald wavevector is consistent.
function dipole_dipole_moment_system(sys::System)
    if is_entangled(sys)
        (; uncontracted, unit_map) = get_entanglement(sys)
        return (uncontracted, unit_map.atom_to_unit)
    else
        return (sys, [UnitPart(i, 1) for i in 1:natoms(sys.crystal)])
    end
end

# Compute Stevens coefficients in the local reference frame
function swt_data(sys::System{0}, measure)
    Na = nsites(sys)
    Nobs = num_observables(measure)

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
    obs = reshape(measure.operators, 1, Nobs, Na)
    obs_localized = [Rs[i]' * obs[1, μ, i] for μ in 1:Nobs, i in 1:Na]

    # Precompute transformed exchange matrices and store in sys.interactions_union.
    for (i, int) in enumerate(sys.interactions_union)
        for c in eachindex(int.pair)
            (; bond, scalar, bilin, biquad, general) = int.pair[c]

            @assert i == bond.i
            j = bond.j

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

            int.pair[c] = PairCoupling(bond, scalar, bilin, biquad, general)
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

# i is a site index for the flattened swt.sys. `measure` was originally
# constructed for a system with nontrivial lattice dims. Use fld1(i, prod(dims))
# to get an atom index for one cell of the unflattened sys.
function formfactor_for_flattened_sys(measure, μ, i)
    sys_dims = size(measure.operators)[3:5]
    measure.formfactors[1, μ, fld1(i, prod(sys_dims))]
end

function observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)
    r_reshaped = sys.crystal.positions[i]
    ff = formfactor_for_flattened_sys(measure, μ, i)
    cis(2π * dot(q_reshaped, r_reshaped)) * compute_form_factor(ff, norm2(q_global))
end

# Linearize each Fourier observable Â(q) = Σᵢ exp(+i q⋅rᵢ) Âᵢ as a coefficient
# vector u in the Nambu basis of Holstein-Primakoff bosons [a_q, a_-q†].
#
# For standard SWT each site has one contribution (ncontribs=1) with zero
# offset. For entangled-unit SWT, SWTDataSUN is populated with ncontribs =
# atoms_per_unit parts that carry distinct position offsets (relative to the
# unit center) and distinct form factor atom indices.
function set_swt_observable_vectors!(u, swt::SpinWaveTheory, q_reshaped, q_global)
    (; sys, measure, data) = swt
    Na = nsites(sys)
    Nobs = num_observables(measure)
    Nf = nflavors(swt)
    L = Nf * Na

    if sys.mode == :SUN
        (; obs_parts, observable_buf) = data::SWTDataSUN
        N = sys.Ns[1]
        natoms_orig = size(measure.operators, 6)
        ncells = div(Na, natoms_orig)
        for μ in 1:Nobs, i in 1:Na
            atom_idx = fld1(i, ncells)
            r = sys.crystal.positions[i]
            fill!(observable_buf, 0)
            for k in axes(obs_parts, 1)
                offset = measure.offsets[k, atom_idx]
                ff = measure.formfactors[k, μ, atom_idx]
                pref = cis(2π * dot(q_reshaped, r + offset)) *
                       compute_form_factor(ff, norm2(q_global))
                observable_buf .+= pref .* obs_parts[k, μ, i]
            end
            for f in 1:Nf
                u[f + (i-1)*Nf,     μ] = observable_buf[f, N]
                u[f + (i-1)*Nf + L, μ] = observable_buf[N, f]
            end
        end
    else
        @assert sys.mode in (:dipole, :dipole_uncorrected)
        (; sqrtS, observables_localized) = data::SWTDataDipole
        for μ in 1:Nobs, i in 1:Na
            pref = observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)
            O = observables_localized[μ, i]
            u[i, μ]   = pref * (sqrtS[i] / √2) * (O[1] + im*O[2])
            u[i+L, μ] = pref * (sqrtS[i] / √2) * (O[1] - im*O[2])
        end
    end

    return u
end
