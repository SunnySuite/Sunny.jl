struct SWTDataDipole
    local_rotations       :: Vector{Mat3}             # Rotations from global to quantization frame
    observables           :: Array{Vec3, 2}           # Observables rotated to local frame (nobs × nsites)
    stevens_coefs         :: Vector{StevensExpansion} # Rotated onsite coupling as Steven expansion
    sqrtS                 :: Vector{Float64}          # Square root of spin magnitudes
    interactions          :: Vector{Interactions}     # Pair couplings rotated to local frames
end


struct SWTDataSUN
    local_unitaries  :: Vector{Matrix{ComplexF64}}    # Transformations from global to quantization frame
    observables      :: Array{HermitianC64, 3}        # Rotated observables (nobs × nunits × nparts)
    spin_ops         :: Array{HermitianC64, 2}        # Spin dipoles in local frame (3 × nbareatoms)
    interactions     :: Vector{Interactions}          # Couplings in local frames, Zeeman in onsite
end

# To facilitate sharing some code with SpinWaveTheorySpiral
abstract type AbstractSpinWaveTheory end

# A correction to the quadratic Hamiltonian, as boson monomials, which defines a
# Gaussian state other than the harmonic vacuum. A correction that was solved
# self-consistently records its `grid`, on which every 1/s correction expanded
# about it must be evaluated, and the outcome of the iteration in `info`.
struct GaussianState
    correction :: Vector{BosonMonomial{2}}
    grid       :: Union{Nothing, BZGrid}
    info       :: NamedTuple
end

GaussianState() = GaussianState(BosonMonomial{2}[], nothing, (;))

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

A `SpinWaveTheory` may carry a Gaussian state other than the harmonic vacuum,
e.g. one returned by [`Sunny.self_consistent`](@ref). Its quadratic Hamiltonian
then includes the state's correction, which determines `dispersion`,
`intensities`, and every other LSWT quantity, while the 1/s corrections treat
the state as the vacuum about which the bosons are expanded.
"""
struct SpinWaveTheory <: AbstractSpinWaveTheory
    sys            :: System
    data           :: Union{SWTDataDipole, SWTDataSUN}
    measure        :: MeasureSpec
    regularization :: Float64
    classical_energy :: Float64
    # The Gaussian state about which the bosons are expanded
    state          :: GaussianState
    # The quadratic Hamiltonian of `state`, and of LSWT itself, Ewald excluded
    quadratic_table :: QuadraticTable
    bare_table      :: QuadraticTable
end

function SpinWaveTheory(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8, energy_ϵ=nothing)
    if !isnothing(energy_ϵ)
        @warn "Keyword argument energy_ϵ is deprecated! Use `regularization` instead."
        regularization = energy_ϵ
    end

    measure = @something measure empty_measurespec(sys)
    if num_lattice_dims(measure) != sys.dims || num_units_per_cell(measure) != natoms(sys.crystal)
        error("Size mismatch. Check that measure is built using consistent system.")
    end

    # Create a new system with dims (1,1,1). A clone happens in all cases.
    new_cryst = resize_and_flatten_crystal(sys.crystal, sys.dims)
    sys = reshape_supercell_aux(sys, new_cryst, (1, 1, 1))

    classical_energy = energy_per_site(sys)

    # Rotate local operators to quantization axis
    data = swt_data(sys, measure)

    L = (sys.mode == :SUN ? sys.Ns[1] - 1 : 1) * natoms(sys.crystal)
    table = QuadraticTable(quadratic_monomials(sys, data), L)
    return SpinWaveTheory(sys, data, measure, regularization, classical_energy, GaussianState(), table, table)
end

# The same theory about the Gaussian state `state`. Everything but the quadratic
# table is shared.
function SpinWaveTheory(swt::SpinWaveTheory, state::GaussianState)
    (; sys, data, measure, regularization, classical_energy, bare_table) = swt
    table = isempty(state.correction) ? bare_table :
        QuadraticTable([quadratic_monomials(sys, data); state.correction], nbands(swt))
    return SpinWaveTheory(sys, data, measure, regularization, classical_energy, state, table, bare_table)
end


function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    printstyled(io, "SpinWaveTheory ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", natoms(swt.sys.crystal), " atoms")
    (; correction, grid, info) = swt.state
    if !isempty(correction)
        println(io, "  Gaussian state: ", isempty(info) ? "user correction" : join(("$k = $v" for (k, v) in pairs(info)), ", "),
                isnothing(grid) ? "" : ", grid $(grid.dims)")
    end
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

# The quadratic table, then the Ewald term, which is dense in q and so not
# compiled, then a small constant shift for positive-definiteness
function dynamical_matrix!(H, swt::SpinWaveTheory, q_reshaped)
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0
    accum_lower!(H, swt.quadratic_table, q_reshaped)
    hermitian_from_lower!(H)

    if !isnothing(uncontracted_system(swt.sys).ewald)
        swt.sys.mode == :SUN ? accum_ewald_SUN!(H, swt, q_reshaped) : accum_ewald_dipole!(H, swt, q_reshaped)
        # Hermitian up to round-off errors, then exactly
        @assert maxdiff(H, H') < 1e-12
        hermitianpart!(H)
    end

    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

function mul_dynamical_matrix(swt, x, qs_reshaped)
    y = zero(x)
    mul_dynamical_matrix!(swt, y, x, qs_reshaped)
end

# Matrix-free y = H x at many wavevectors, where row q of x and y belongs to
# `qs_reshaped[q]`
function mul_dynamical_matrix!(swt, y, x, qs_reshaped; phases=zeros(ComplexF64, length(qs_reshaped)))
    @assert size(x) == size(y) == (length(qs_reshaped), 2nbands(swt))
    isnothing(swt.sys.ewald) || error("Ewald not supported")
    y .= 0
    accum_mul!(y, x, swt.quadratic_table, qs_reshaped, phases)
    @. y += swt.regularization * x
    return y
end


# Prepare local operators, observables, and couplings for spin wave calculation
# by rotating into the local reference frame as defined by the ground state. The
# couplings of `sys` itself are left in the global frame.
function swt_data(sys::System{N}, @nospecialize measure) where N
    # Calculate transformation matrices into local reference frames
    Na = nsites(sys)
    Nb = nbaresites(sys)
    nparts = num_parts_per_unit(measure)
    Nobs = num_observables(measure)
    flat_ops = reshape(measure.observables::Array{HermitianC64, 6}, Nobs, Na, nparts)

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = Vector{Matrix{ComplexF64}}(undef, Na)
    observables = Array{HermitianC64}(undef, Nobs, Na, nparts)
    spin_ops = fill(Hermitian(zeros(ComplexF64, N, N)), 3, Nb)
    onsites = [int.onsite for int in sys.interactions_union]

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
        for p in 1:nparts, μ in 1:Nobs
            observables[μ, i, p] = Hermitian(U' * flat_ops[μ, i, p] * U)
        end

        for ai in atoms_in_unit(sys, i)
            # Incorporate Zeeman coupling into onsite
            S = spin_ops_embedded(sys, ai)
            (; gs, extfield) = uncontracted_system(sys)
            onsites[i] += (extfield[ai]' * gs[ai]) * S

            # Store spin dipoles in the local frame
            for α in 1:3
                spin_ops[α, ai] = Hermitian(U' * S[α] * U)
            end
        end
    end

    interactions = map(1:Na) do i
        Ui = local_unitaries[i]
        int = sys.interactions_union[i]

        # Rotate onsite coupling
        onsite = Hermitian(Ui' * onsites[i] * Ui)

        # Convert each pair coupling to a compressed tensor decomposed form and
        # rotate it into the local reference frames of both sites.
        pair_new = PairCoupling[]
        for pc in int.pair
            (; bond) = pc
            @assert bond.i == i

            N1, N2 = sys.Ns[bond.i], sys.Ns[bond.j]
            Uj = local_unitaries[bond.j]
            data_rotated = map(pair_coupling_tensor_data(pc, N1, N2)) do (A, B)
                (Hermitian(Ui'*A*Ui), Hermitian(Uj'*B*Uj))
            end

            gen1 = spin_matrices_of_dim(; N=N1)
            gen2 = spin_matrices_of_dim(; N=N2)
            general_rotated = TensorDecomposition(gen1, gen2, data_rotated)

            push!(pair_new, PairCoupling(bond, 0.0, 0.0, 0.0, general_rotated))
        end
        Interactions(onsite, pair_new)
    end

    return SWTDataSUN(local_unitaries, observables, spin_ops, interactions)
end


function swt_data(sys::System{0}, @nospecialize measure)
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
    obs = reshape(measure.observables::Array{Vec3, 6}, Nobs, Na)
    obs_localized = [Rs[i]' * obs[μ, i] for μ in 1:Nobs, i in 1:Na]

    # Precompute transformed exchange matrices. The onsite coupling is kept in
    # the global frame, its rotated form being `stevens_coefs`.
    interactions = map(sys.interactions_union) do int
        pair = map(int.pair) do pc
            (; bond, scalar, bilin, biquad, general) = pc
            i = bond.i

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

            PairCoupling(bond, scalar, bilin, biquad, general)
        end
        Interactions(int.onsite, pair)
    end

    # Rotated Stevens expansion
    cs = map(sys.interactions_union, Rs) do int, R
        iszero(R) ? empty_anisotropy(sys.mode, 0) : rotate_operator(int.onsite, R)
    end

    # Square root of spin magnitudes
    sqrtS = sqrt.(vec(sys.κs))

    return SWTDataDipole(Rs, obs_localized, cs, sqrtS, interactions)
end

# Fourier prefactor for observable μ acting on a flattened unit i of swt.sys,
# with optional part index.
function observable_prefactor(measure, μ, i, q_reshaped, q_global, sys; part=1)
    unit = fld1(i, prod(num_lattice_dims(measure)))
    r = sys.crystal.positions[i] + measure.offsets[unit, part] ./ Vec3(num_lattice_dims(measure))
    ff = measure.formfactors[μ, unit, part]
    cis(2π * dot(q_reshaped, r)) * compute_form_factor(ff, norm2(q_global))
end

# Linearize each Fourier observable Â(q) = Σᵢ exp(+i q⋅rᵢ) Âᵢ as a coefficient
# vector u in the Nambu basis of Holstein-Primakoff bosons [a_q, a_-q†].
function set_swt_observable_vectors!(u, swt::SpinWaveTheory, q_reshaped, q_global)
    (; sys, measure, data) = swt
    Na = nsites(sys)
    Nobs = num_observables(measure)
    Nf = nflavors(swt)
    L = Nf * Na

    if sys.mode == :SUN
        (; observables) = data::SWTDataSUN
        @assert allequal(sys.Ns)
        N = first(sys.Ns)
        for μ in 1:Nobs, i in 1:Na
            view(u, (i-1)*Nf .+ (1:Nf), μ) .= 0
            view(u, (i-1)*Nf .+ (1:Nf) .+ L, μ) .= 0
            for p in 1:num_parts_per_unit(measure)
                pref = observable_prefactor(measure, μ, i, q_reshaped, q_global, sys; part=p)
                O = observables[μ, i, p]
                for f in 1:Nf
                    u[f + (i-1)*Nf,     μ] += pref * O[f, N]
                    u[f + (i-1)*Nf + L, μ] += pref * O[N, f]
                end
            end
        end
    else
        @assert sys.mode in (:dipole, :dipole_uncorrected)
        num_parts_per_unit(measure) == 1 || error("Distinct parts unsupported in dipole-mode.")
        (; sqrtS, observables) = data::SWTDataDipole
        for μ in 1:Nobs, i in 1:Na
            pref = observable_prefactor(measure, μ, i, q_reshaped, q_global, sys)
            O = observables[μ, i]
            u[i, μ]   = pref * (sqrtS[i] / √2) * (O[1] + im*O[2])
            u[i+L, μ] = pref * (sqrtS[i] / √2) * (O[1] - im*O[2])
        end
    end

    return u
end
