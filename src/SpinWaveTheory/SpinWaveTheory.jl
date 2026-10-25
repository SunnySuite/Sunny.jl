# A spin wave theory that reports `intensities_bands`, from which the broadened
# and static intensities follow. Implemented by SpinWaveTheory and
# SpinWaveTheorySpiral.
abstract type AbstractSpinWaveTheory end

struct SWTDataDipole
    local_rotations       :: Vector{Mat3}             # Rotations from global to quantization frame
    # Observables rotated to local frame (nobs × nsites), each a spin component
    # 𝐎⋅𝐒 or, for a measure of general observables, a Stevens expansion
    observables           :: Union{Array{Vec3, 2}, Array{StevensExpansion, 2}}
    stevens_coefs         :: Vector{StevensExpansion} # Rotated onsite coupling as Steven expansion
    sqrtS                 :: Vector{Float64}          # Square root of spin magnitudes
    interactions          :: Vector{Interactions}     # Pair couplings rotated to local frames
end


struct SWTDataSUN
    local_unitaries  :: Vector{Matrix{ComplexF64}}    # Transformations from global to quantization frame
    observables      :: Array{HermitianC64, 3}        # Rotated observables (nobs × nsites × nparts)
    spin_ops         :: Array{HermitianC64, 2}        # Spin dipoles in local frame (3 × bare sites)
    interactions     :: Vector{Interactions}          # Couplings in local frames, Zeeman in onsite
end

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

# The monomials of H₃ and H₄, built on first use. They depend on the structure
# and couplings but not on the Gaussian state, so every state of one theory
# shares them. Only the 1/s corrections need them, and building them costs
# several times more than the rest of `SpinWaveTheory`. Consumers must not
# mutate the lists returned.
mutable struct MonomialCache
    const lock :: ReentrantLock
    cubic      :: Union{Nothing, Vector{BosonMonomial{3}}}
    quartic    :: Union{Nothing, Vector{BosonMonomial{4}}}
end

MonomialCache() = MonomialCache(ReentrantLock(), nothing, nothing)

# The field `name` of `cache`, built by `f` if it is not yet present
function cached!(f, cache::MonomialCache, name::Symbol)
    @lock cache.lock begin
        v = getfield(cache, name)
        isnothing(v) ? setfield!(cache, name, f()) : v
    end
end

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
e.g. one returned by `Sunny.self_consistent`. Its quadratic Hamiltonian
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
    monomials       :: MonomialCache
    # One-boson words of the observables; see `observable_coefficients`
    observable_coefs :: Array{ComplexF64, 3}
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

    # A clone protects against later mutation of `sys`. The quasi-particle
    # bands are indexed by the sites of the whole system, and the wavevectors
    # by the reciprocal lattice of its magnetic cell; see `to_supercell_rlu`.
    sys = clone_system(sys)
    classical_energy = energy_per_site(sys)

    # Rotate local operators to quantization axis
    data = swt_data(sys, measure)

    return SpinWaveTheory(sys, data, measure, regularization, classical_energy)
end

# The harmonic theory of the structure in `sys`, given its local `data`
function SpinWaveTheory(sys::System, data, measure, regularization, classical_energy)
    table = QuadraticTable(quadratic_monomials(sys, data), nbands(sys))
    coefs = observable_coefficients(sys, data, measure)
    return SpinWaveTheory(sys, data, measure, regularization, classical_energy, GaussianState(), table, table,
                          MonomialCache(), coefs)
end

# The same theory about the Gaussian state `state`. Everything but the quadratic
# table is shared.
function SpinWaveTheory(swt::SpinWaveTheory, state::GaussianState)
    (; sys, data, measure, regularization, classical_energy, bare_table, monomials, observable_coefs) = swt
    table = isempty(state.correction) ? bare_table :
        QuadraticTable([quadratic_monomials(sys, data); state.correction], nbands(swt))
    return SpinWaveTheory(sys, data, measure, regularization, classical_energy, state, table, bare_table,
                          monomials, observable_coefs)
end


function Base.show(io::IO, ::MIME"text/plain", swt::SpinWaveTheory)
    printstyled(io, "SpinWaveTheory ", mode_to_str(swt.sys), "\n"; bold=true, color=:underline)
    println(io, "  ", nsites(swt.sys), " sites")
    (; correction, grid, info) = swt.state
    if !isempty(correction)
        println(io, "  Gaussian state: ", isempty(info) ? "user correction" : join(("$k = $v" for (k, v) in pairs(info)), ", "),
                isnothing(grid) ? "" : ", grid $(grid.dims)")
    end
end

nflavors(sys::System) = sys.mode == :SUN ? sys.Ns[1]-1 : 1
nbands(sys::System) = nflavors(sys) * nsites(sys)
nflavors(swt::SpinWaveTheory) = nflavors(swt.sys)
nbands(swt::SpinWaveTheory) = nbands(swt.sys)

# For a site with linear index i in `1:nsites(sys)`, return the linear index j
# of the site that `bond` reaches, and the number of magnetic cells by which it
# wraps.
function bond_target(sys::System, i::Int, bond::Bond)
    cell = Tuple(to_cell(eachsite(sys)[i])) .- 1 .+ Tuple(bond.n)
    j = LinearIndices(sys.dipoles)[(mod.(cell, sys.dims) .+ 1)..., bond.j]
    return (j, Vec3(fld.(cell, sys.dims)))
end

# The interactions of the site with linear index i in `1:nsites(sys)`, stored
# per sublattice in a homogeneous system
interactions_at(sys::System, i::Int) =
    is_homogeneous(sys) ? sys.interactions_union[to_atom(eachsite(sys)[i])] : sys.interactions_union[i]

# Wavevectors in RLU of the original crystal, `q`, and of the magnetic cell of
# `sys`, `q_supercell`. The magnetic cell spans `sys.dims` cells of
# `sys.crystal`, which may itself be a reshaped cell.
function from_supercell_rlu(sys::System, q_supercell)
    return orig_crystal(sys).recipvecs \ (sys.crystal.recipvecs * (q_supercell ./ Vec3(sys.dims)))
end

function to_supercell_rlu(sys::System, q)
    return Vec3(sys.dims) .* (sys.crystal.recipvecs \ (orig_crystal(sys).recipvecs * q))
end

function dynamical_matrix(swt::SpinWaveTheory, q_supercell)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    dynamical_matrix!(H, swt, q_supercell)
    return Hermitian(H)
end

# The quadratic table, then the Ewald term, which is dense in q and so not
# compiled, then a small constant shift for positive-definiteness
function dynamical_matrix!(H, swt::SpinWaveTheory, q_supercell)
    L = nbands(swt)
    @assert size(H) == (2L, 2L)
    H .= 0
    accum_lower!(H, swt.quadratic_table, q_supercell)
    hermitian_from_lower!(H)

    if !isnothing(uncontracted_system(swt.sys).ewald)
        accum_ewald!(H, swt, q_supercell)
        # Hermitian up to round-off errors, then exactly
        @assert maxdiff(H, H') < 1e-12
        hermitianpart!(H)
    end

    for i in 1:2L
        H[i, i] += swt.regularization
    end
end

function mul_dynamical_matrix(swt, x, qs_supercell)
    y = zero(x)
    mul_dynamical_matrix!(swt, y, x, qs_supercell)
end

# Matrix-free y = H x at many wavevectors, where row q of x and y belongs to
# `qs_supercell[q]`
function mul_dynamical_matrix!(swt, y, x, qs_supercell; phases=zeros(ComplexF64, length(qs_supercell)))
    @assert size(x) == size(y) == (length(qs_supercell), 2nbands(swt))
    isnothing(uncontracted_system(swt.sys).ewald) || error("Ewald not supported")
    y .= 0
    accum_mul!(y, x, swt.quadratic_table, qs_supercell, phases)
    @. y += swt.regularization * x
    return y
end


# Prepare local operators, observables, and couplings for spin wave calculation
# by rotating into the local reference frame as defined by the ground state. The
# couplings of `sys` itself are left in the global frame.
function swt_data(sys::System{N}, @nospecialize measure) where N
    # Calculate transformation matrices into local reference frames
    Na = nsites(sys)
    nparts = num_parts_per_unit(measure)
    Nobs = num_observables(measure)
    flat_ops = reshape(measure.observables::Array{HermitianC64, 6}, Nobs, Na, nparts)

    # Preallocate buffers for local unitaries and observables. Spin dipoles are
    # stored per bare site of the uncontracted system.
    usys = uncontracted_system(sys)
    local_unitaries = Vector{Matrix{ComplexF64}}(undef, Na)
    observables = Array{HermitianC64}(undef, Nobs, Na, nparts)
    spin_ops = fill(Hermitian(zeros(ComplexF64, N, N)), 3, nsites(usys))

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

        # Store spin dipoles of each bare site in the local frame
        for (bi, _) in bare_indices_and_cells(sys, i)
            S = spin_ops_embedded(sys, to_atom(eachsite(usys)[bi]))
            for α in 1:3
                spin_ops[α, bi] = Hermitian(U' * S[α] * U)
            end
        end
    end

    interactions = map(1:Na) do i
        Ui = local_unitaries[i]
        int = interactions_at(sys, i)

        # Rotate onsite coupling, with the Zeeman coupling of each bare site
        onsite = int.onsite
        for (bi, _) in bare_indices_and_cells(sys, i)
            S = spin_ops_embedded(sys, to_atom(eachsite(usys)[bi]))
            onsite += (usys.extfield[bi]' * usys.gs[bi]) * S
        end
        onsite = Hermitian(Ui' * onsite * Ui)

        # Convert each pair coupling to a compressed tensor decomposed form and
        # rotate it into the local reference frames of both sites.
        pair_new = PairCoupling[]
        for pc in int.pair
            (; bond) = pc
            (j, _) = bond_target(sys, i, bond)

            N1, N2 = sys.Ns[i], sys.Ns[j]
            Uj = local_unitaries[j]
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


# The local `frames` may instead be given, e.g. to continue a family of frames
# smoothly, as `relax_frames` does. Column 3 of each must be the direction of
# its dipole.
function swt_data(sys::System{0}, @nospecialize measure; frames=nothing)
    Na = nsites(sys)
    Nobs = num_observables(measure)

    # Operators for rotating vectors into local frame
    Rs = !isnothing(frames) ? frames : map(1:Na) do i
        # Defines quantization axis.
        n = sys.dipoles[i]

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
    # in practice left-multiply column vector by R'. A Stevens expansion rotates
    # as the onsite coupling does.
    num_parts_per_unit(measure) == 1 || error("Distinct parts unsupported in dipole-mode.")
    obs = reshape(measure.observables, Nobs, Na)
    obs_localized = if obs isa Array{Vec3}
        [Rs[i]' * obs[μ, i] for μ in 1:Nobs, i in 1:Na]
    else
        [iszero(Rs[i]) ? obs[μ, i] * 0.0 : rotate_operator(obs[μ, i], Rs[i]) for μ in 1:Nobs, i in 1:Na]
    end

    # Precompute transformed exchange matrices, per site. The onsite coupling
    # is kept in the global frame, its rotated form being `stevens_coefs`.
    interactions = map(1:Na) do i
        int = interactions_at(sys, i)
        pair = map(int.pair) do pc
            (; bond, scalar, bilin, biquad, general) = pc
            (j, _) = bond_target(sys, i, bond)

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
    cs = map(1:Na, Rs) do i, R
        iszero(R) ? empty_anisotropy(sys.mode, 0) : rotate_operator(interactions_at(sys, i).onsite, R)
    end

    # Square root of spin magnitudes
    sqrtS = sqrt.(vec(sys.κs))

    return SWTDataDipole(Rs, obs_localized, cs, sqrtS, interactions)
end

# Fourier prefactor for observable μ acting on site i of `sys`, with optional
# part index. The measure offsets are fractional coordinates of the cell.
function observable_prefactor(measure, μ, i, q_global, sys; part=1)
    site = eachsite(sys)[i]
    r = global_position_at(sys, site) + sys.crystal.latvecs * measure.offsets[to_atom(site), part]
    ff = measure.formfactors[μ, to_atom(site), part]
    cis(dot(q_global, r)) * compute_form_factor(ff, norm2(q_global))
end

# Coefficients C[a, μ, p] of the one-boson words of part p of each local
# observable μ, by Nambu label a, in the local frame and without the Fourier
# phase or form factor that `observable_prefactor` supplies. These are the
# amplitudes of LSWT, the leading order of the words of
# `observable_monomials`.
function observable_coefficients(sys::System, data, measure::MeasureSpec)
    C = zeros(ComplexF64, 2nbands(sys), num_observables(measure), num_parts_per_unit(measure))
    for p in axes(C, 3), μ in axes(C, 2), i in 1:nsites(sys)
        for (; c, as) in observable_monomials(sys, data, μ, i, Val{1}(); part=p, beyond_tree=false)
            C[as[1], μ, p] += c
        end
    end
    return C
end

# Linearize each Fourier observable Â(q) = Σᵢ exp(+i q⋅rᵢ) Âᵢ as a coefficient
# vector u in the Nambu basis of Holstein-Primakoff bosons [a_q, a_-q†].
#
# The fields `sys` and `measure` of `SpinWaveTheory` are abstract, to limit
# specialization, so they are passed through a function barrier. Otherwise each
# observable would dispatch dynamically.
observable_vectors!(u, swt::SpinWaveTheory, q_global) =
    observable_vectors!(u, swt, swt.sys, swt.measure, q_global)

function observable_vectors!(u, swt::SpinWaveTheory, sys::System, measure::MeasureSpec, q_global)
    fill!(u, 0)
    return accum_observable_vectors!(u, sys, measure, swt.observable_coefs, q_global)
end

# Adds to `u` the observables whose one-boson words have coefficients C[a, μ, p],
# as in `observable_coefficients`. The coefficient of the operator labeled a is
# stored in u[ā, μ], where ā = mod1(a+L, 2L), and each part p of a unit carries
# its own Fourier phase and form factor.
function accum_observable_vectors!(u, sys::System, measure::MeasureSpec, C, q_global)
    L = nbands(sys)
    Nf = nflavors(sys)
    for p in axes(C, 3), μ in axes(C, 2), i in 1:nsites(sys)
        pref = observable_prefactor(measure, μ, i, q_global, sys; part=p)
        for a in (i-1)*Nf .+ (1:Nf)
            u[a, μ]   += pref * C[L+a, μ, p]
            u[L+a, μ] += pref * C[a, μ, p]
        end
    end
    return u
end

# Returns |1 + nB(ω)| where nB(ω) = 1 / (exp(βω) - 1) is the Bose function.
# Equivalent to |1 / expm1(-βω)| where expm1(x) = e^x-1.
function thermal_prefactor(ω; kT)
    @assert kT >= 0
    iszero(ω) && return Inf
    return abs(1 / expm1(-ω/kT))
end

"""
    excitations!(T, tmp, swt::SpinWaveTheory, q)

Given a wavevector `q`, solves for the matrix `T` representing quasi-particle
excitations, and returns a list of quasi-particle energies. Both `T` and `tmp`
must be supplied as ``2L×2L`` complex matrices, where ``L`` is the number of
bands for a single ``𝐪`` value.

The columns of `T` are understood to be contracted with the Holstein-Primakoff
bosons ``[𝐛_𝐪, 𝐛_{-𝐪}^†]``. The first ``L`` columns provide the eigenvectors
of the quadratic Hamiltonian for the wavevector ``𝐪``. The next ``L`` columns
of `T` describe eigenvectors for ``-𝐪``. The return value is a vector with
similar grouping: the first ``L`` values are energies for ``𝐪``, and the next
``L`` values are the _negation_ of energies for ``-𝐪``.

    excitations!(T, tmp, swt::SpinWaveTheorySpiral, q; branch)

Calculations on a [`SpinWaveTheorySpiral`](@ref) additionally require a `branch`
index. The possible branches ``(1, 2, 3)`` correspond to scattering processes
``𝐪 - 𝐤, 𝐪, 𝐪 + 𝐤`` respectively, where ``𝐤`` is the ordering wavevector.
Each branch will contribute ``L`` excitations, where ``L`` is the number of
spins in the magnetic cell. This yields a total of ``3L`` excitations for a
given momentum transfer ``𝐪``.
"""
function excitations!(T, tmp, swt::SpinWaveTheory, q)
    L = nbands(swt)
    size(T) == size(tmp) == (2L, 2L) || error("Arguments T and tmp must be $(2L)×$(2L) matrices")
    ws = BogoliubovWorkspace(L)
    energies = excitations!(ws, tmp, swt, q)
    copyto!(T, ws.T)
    return energies
end

# Allocation-free variant that stores the transformation in `ws.T`. The returned
# energies alias `ws.ε`.
excitations!(ws::BogoliubovWorkspace, H, swt::SpinWaveTheory, q) =
    excitations_supercell!(ws, H, swt, to_supercell_rlu(swt.sys, q))

# The same at a wavevector in supercell RLU. A quadratic Hamiltonian that is not
# positive definite has no vacuum. For LSWT that means the structure is not an
# energy minimum; for a corrected Gaussian state, the correction is to blame.
function excitations_supercell!(ws::BogoliubovWorkspace, H, swt::SpinWaveTheory, q_supercell)
    dynamical_matrix!(H, swt, q_supercell)
    try
        return bogoliubov!(ws, H)
    catch err
        err isa PosDefException || rethrow()
        q = vec3_to_string(from_supercell_rlu(swt.sys, q_supercell))
        rethrow(InstabilityError(isempty(swt.state.correction) ?
            "Not an energy-minimum; wavevector q = $q unstable." :
            "Quadratic Hamiltonian of the Gaussian state is not positive definite at wavevector q = $q. If \
             quantum fluctuations stabilize the structure, see `Sunny.self_consistent`."))
    end
end

"""
    excitations(swt::SpinWaveTheory, q)
    excitations(swt::SpinWaveTheorySpiral, q; branch)

Returns a pair `(energies, T)` providing the excitation energies and
eigenvectors. Prefer [`excitations!`](@ref) for performance, which avoids matrix
allocations. See the documentation of [`excitations!`](@ref) for more details.
"""
function excitations(swt::SpinWaveTheory, q)
    L = nbands(swt)
    T = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    energies = excitations!(T, copy(H), swt, q)
    return (energies, T)
end

"""
    dispersion(swt::SpinWaveTheory, qpts)

Given a list of wavevectors `qpts` in reciprocal lattice units (RLU), returns
excitation energies for each band. The return value `ret` is 2D array, and
should be indexed as `ret[band_index, q_index]`.
"""
function dispersion(swt::SpinWaveTheory, qpts)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    disp = zeros(L, length(qpts.qs))
    for (iq, q) in enumerate(qpts.qs)
        view(disp, :, iq) .= view(excitations(swt, q)[1], 1:L)
    end
    return reshape(disp, L, size(qpts.qs)...)
end

"""
    intensities_bands(swt::SpinWaveTheory, qpts; kT=0, threaded=false)

Calculate spin wave excitation bands for a set of ``𝐪``-points in reciprocal
space. This calculation is analogous to [`intensities`](@ref), but does not
perform line broadening of the bands. Use [`load_blas_for_threading`](@ref) and
set `threaded=true` to parallelize the calculation using Julia threads.
"""
function intensities_bands(swt::SpinWaveTheory, qpts; kT=0, with_negative=false, threaded=false)
    num_observables(swt.measure) == 0 && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    with_negative && error("Option `with_negative=true` not yet supported.")
    # The fields `sys` and `measure` are abstract, to limit specialization, so
    # the loop over wavevectors is behind a function barrier
    return intensities_bands(swt, swt.sys, swt.measure, convert(AbstractQPoints, qpts); kT, threaded)
end

function intensities_bands(swt::SpinWaveTheory, sys::System, measure::MeasureSpec, qpts::AbstractQPoints; kT, threaded)
    cryst = orig_crystal(sys)

    # Number of (magnetic) atoms in magnetic cell
    Na = nsites(sys)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Temporary storage for pair correlations
    Nobs = num_observables(measure)
    Ncorr = num_correlations(measure)

    disp = zeros(Float64, L, Nq)
    data = zeros(eltype(measure), L, Nq)

    newbuf() = (
        ws   = BogoliubovWorkspace(L),
        H    = zeros(ComplexF64, 2L, 2L),
        u    = zeros(ComplexF64, 2L, Nobs),
        Avec = zeros(ComplexF64, Nobs),
        corr = zeros(ComplexF64, Ncorr),
    )

    function calc_iq!(buf, iq)
        (; ws, H, u, Avec, corr) = buf
        T = ws.T
        q = qpts.qs[iq]
        q_supercell = to_supercell_rlu(sys, q)
        q_global = cryst.recipvecs * q
        view(disp, :, iq) .= view(excitations_supercell!(ws, H, swt, q_supercell), 1:L)

        # Linearized observables Â_μ(q) in Holstein-Primakoff bosons
        observable_vectors!(u, swt, sys, measure, q_global)

        # Fill `intensity` array
        for band in 1:L
            for μ in 1:Nobs
                # Left matrix amplitude ⟨0|Â_μ(q)†|n⟩
                Avec[μ] = dot(view(u, :, μ), view(T, :, band))
            end
            map!(corr, measure.corr_pairs) do (μ, ν)
                # Pair correlations ⟨0|Â_μ(q)†|n⟩⟨n|Â_ν(q)|0⟩
                Avec[μ] * conj(Avec[ν]) / Ncells
            end
            data[band, iq] = thermal_prefactor(disp[band, iq]; kT) * measure.combiner(q_global, corr)
        end
    end

    foreach_chunked(calc_iq!, newbuf, 1:Nq; threaded, warn_blas=true)

    # Reassigning `disp` and `data` would box these captured variables
    disp_reshaped = reshape(disp, L, size(qpts.qs)...)
    data_reshaped = reshape(data, L, size(qpts.qs)...)

    return BandIntensities(cryst, qpts, disp_reshaped, data_reshaped)
end

"""
    intensities!(data, swt::SpinWaveTheory, qpts; energies, kernel, kT=0, threaded=false)
    intensities!(data, sc::SampledCorrelations, qpts; energies, kernel=nothing, kT=0)

Like [`intensities`](@ref), but makes use of storage space `data` to avoid
allocation costs.
"""
function intensities!(data, swt::AbstractSpinWaveTheory, qpts; energies, kernel::AbstractBroadening, kT=0, threaded=false)
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    bands = intensities_bands(swt, qpts; kT, threaded)
    @assert eltype(bands) == eltype(data)
    broaden!(data, bands; energies, kernel, threaded)
    return Intensities(bands.crystal, bands.qpts, collect(Float64, energies), data)
end

"""
    intensities(swt::SpinWaveTheory, qpts; energies, kernel, kT=0, threaded=false)
    intensities(sc::SampledCorrelations, qpts; energies, kernel=nothing, kT)

Calculates dynamical pair correlation intensities for a set of ``𝐪``-points in
reciprocal space.

Linear spin wave theory calculations are performed with an instance of
[`SpinWaveTheory`](@ref). The alternative [`SpinWaveTheorySpiral`](@ref) allows
to study generalized spiral orders with a single, incommensurate-``𝐤`` ordering
wavevector. Another alternative [`SpinWaveTheoryKPM`](@ref) is favorable for
calculations on large magnetic cells, and allows to study systems with disorder.
An optional nonzero temperature `kT` will scale intensities by the quantum
thermal occupation factor ``|1 + n_B(ω)|`` where ``n_B(ω) = 1/(e^{βω}-1)`` is
the Bose function.

Intensities can also be calculated for `SampledCorrelations` associated with
classical spin dynamics. In this case, thermal broadening will already be
present, and the line-broadening `kernel` may be omitted. Conversely, the
parameter `kT` becomes required. If positive, it will introduce an intensity
correction factor ``|βω (1 + n_B(ω))|`` that undoes the occupation factor for
the classical Boltzmann distribution and applies the quantum thermal occupation
factor. The special choice `kT = nothing` will suppress the classical-to-quantum
correction factor, and yield statistics consistent with the classical Boltzmann
distribution. If a `kernel` is provided, it will be used to perform a
convolution along the energy axis on top of any intrinsic broadening already
present in the correlation data. In this case, `energies` must be an explicit
list specifying the desired output energies, exactly as for `SpinWaveTheory`.
"""
function intensities(swt::AbstractSpinWaveTheory, qpts; energies, kernel::AbstractBroadening, kT=0, threaded=false)
    return broaden(intensities_bands(swt, qpts; kT, threaded); energies, kernel, threaded)
end

"""
    intensities_static(swt::SpinWaveTheory, qpts; bounds=(-Inf, Inf), kernel=nothing, kT=0, threaded=false)
    intensities_static(sc::SampledCorrelations, qpts; bounds=(-Inf, Inf), kT)
    intensities_static(sc::SampledCorrelationsStatic, qpts)

Like [`intensities`](@ref), but integrates the dynamical correlations
``\\mathcal{S}(𝐪, ω)`` over a range of energies ``ω``. By default, the
integration `bounds` are ``(-∞, ∞)``, yielding the instantaneous (equal-time)
correlations.

In [`SpinWaveTheory`](@ref), the integral will be realized as a sum over
discrete bands. Alternative calculation methods are
[`SpinWaveTheorySpiral`](@ref) and [`SpinWaveTheoryKPM`](@ref).

Classical dynamics data in [`SampledCorrelations`](@ref) can also be used to
calculate static intensities. In this case, the domain of integration will be a
finite grid of available `energies`. Here, the parameter `kT` will be used to
account for the quantum thermal occupation of excitations, as documented in
[`intensities`](@ref).

Static intensities calculated from [`SampledCorrelationsStatic`](@ref) are
dynamics-independent. Instead, instantaneous correlations sampled from the
classical Boltzmann distribution will be reported.
"""
function intensities_static(swt::AbstractSpinWaveTheory, qpts; bounds=(-Inf, Inf), kernel=nothing, kT=0, threaded=false)
    res = intensities_bands(swt, qpts; kT, threaded)  # TODO: with_negative=true
    data_reduced = zeros(eltype(res.data), size(res.data)[2:end])
    for iq in CartesianIndices(data_reduced), ib in axes(res.data, 1)
        ϵ = res.disp[ib, iq]
        if isnothing(kernel) || bounds == (-Inf, Inf)
            if bounds[1] <= ϵ < bounds[2]
                data_reduced[iq] += res.data[ib, iq]
            end
        else
            isnothing(kernel.integral) && error("Kernel must provide integral")
            ihi = kernel.integral(bounds[2] - ϵ)
            ilo = kernel.integral(bounds[1] - ϵ)
            data_reduced[iq] += res.data[ib, iq] * (ihi - ilo)
        end
    end
    StaticIntensities(res.crystal, res.qpts, data_reduced)
end
