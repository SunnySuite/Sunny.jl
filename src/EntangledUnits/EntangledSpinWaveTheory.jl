struct SWTDataEntangled
    local_unitaries           :: Vector{Matrix{ComplexF64}}   # Aligns to quantization axis on each site
    observables_localized     :: Array{HermitianC64, 3}   # Observables in local frame for each subsite (for intensity calcs)
    observable_buf            :: Array{ComplexF64, 2}   # Buffer for use while constructing boson rep of observables
end

struct EntangledSpinWaveTheory <: AbstractDirectSpinWaveTheory
    sys              :: System
    data             :: SWTDataEntangled
    measure          :: MeasureSpec
    regularization   :: Float64

    crystal_origin   :: Crystal
    contraction_info :: CrystalContractionInfo
    Ns_unit          :: Vector{Vector{Int64}}
end


function SpinWaveTheory(esys::EntangledSystem; measure, regularization=1e-8)
    (; sys_contracted, sys_uncontracted) = esys
    crystal_origin = orig_crystal(sys_uncontracted)
    if !isnothing(sys_contracted.ewald)
        error("EntangledSpinWaveTheory does not yet support long-range dipole-dipole interactions.")
    end

    measure = @something measure empty_measurespec(sys_contracted)
    if nsites(sys_uncontracted) != prod(size(measure.observables)[2:5])
        error("Size mismatch. Check that measure is built using consistent system.")
    end

    # Reshape (flatten) both sys_contracted and sys_uncontracted in corresponding ways
    sys_contracted, sys_uncontracted = map([sys_contracted, sys_uncontracted]) do sys
        new_shape = cell_shape(sys) * diagm(Vec3(sys.dims))
        new_cryst = reshape_crystal(orig_crystal(sys), new_shape)

        # Create a new system with dims (1,1,1). A clone happens in all cases.
        return reshape_supercell_aux(sys, new_cryst, (1,1,1))
    end

    # Construct the new contraction_info (i.e. reconstruct maps between
    # entangled and unentangled systems)
    new_units = units_for_reshaped_system(sys_uncontracted, esys)
    _, new_contraction_info = contract_crystal(sys_uncontracted.crystal, new_units)
    new_Ns_unit = Ns_in_units(sys_uncontracted, new_contraction_info)

    data = swt_data_entangled(sys_contracted, sys_uncontracted, new_Ns_unit, new_contraction_info, measure)

    return EntangledSpinWaveTheory(sys_contracted, data, measure, regularization, crystal_origin, new_contraction_info, new_Ns_unit)
end

function Base.show(io::IO, ::MIME"text/plain", swt::EntangledSpinWaveTheory)
    printstyled(io, "EntangledSpinWaveTheory\n"; bold=true, color=:underline)
    println(io, "Entangled units in magnetic supercell: $(natoms(swt.sys.crystal))")
    # Add observable info?
end

nbands(swt::EntangledSpinWaveTheory) = (swt.sys.Ns[1]-1)  * natoms(swt.sys.crystal)

function to_reshaped_rlu(esys::EntangledSystem, q)
    (; sys_contracted, sys_uncontracted) = esys
    return sys_contracted.crystal.recipvecs \ (orig_crystal(sys_uncontracted).recipvecs * q)
end

function to_reshaped_rlu(eswt::EntangledSpinWaveTheory, q)
    (; sys, crystal_origin) = eswt
    return sys.crystal.recipvecs \ (crystal_origin.recipvecs * q)
end

# obs are observables _given in terms of `sys_original`_
function swt_data_entangled(sys_contracted, sys_uncontracted, Ns_unit, contraction_info, measure)

    # Calculate transformation matrices into local reference frames
    N = sys_contracted.Ns[1] # Assume uniform contraction for now
    nunits = nsites(sys_contracted)
    natoms = nsites(sys_uncontracted)
    nobs = size(measure.observables, 1)
    observables = reshape(measure.observables, nobs, natoms)

    # Check to make sure all "units" are the same -- this can be generalized later.
    atoms_per_unit = [length(info) for info in contraction_info.inverse]
    @assert allequal(atoms_per_unit) "All units must have the same number of interior sites"
    atoms_per_unit = atoms_per_unit[1]

    Ns_contracted = map(Ns -> prod(Ns), Ns_unit)
    @assert allequal(Ns_contracted) "All units must have the same dimension local Hilbert space"
    @assert Ns_contracted[1] == N "Unit dimension inconsistent with system"  # Sanity check. This should never happen.

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = Vector{Matrix{ComplexF64}}(undef, nunits)
    observables_localized = Array{HermitianC64}(undef, atoms_per_unit, nobs, nunits)
    observable_buf = zeros(ComplexF64, N, N)

    # Find local unitaries (for units)
    for i in 1:nunits
        # Create unitary that rotates [0, ..., 0, 1] into ground state direction
        # Z that defines quantization axis
        Z = sys_contracted.coherents[i]
        U = hcat(nullspace(Z'), Z)
        local_unitaries[i] = U

        # Rotate observables into local reference frames
        for μ in 1:num_observables(measure)
            A = observables[μ, i]
            for ui in 1:atoms_per_unit
                A_prod = local_op_to_product_space(A, ui, Ns_unit[i])
                observables_localized[ui, μ, i] = Hermitian(U' * convert(Matrix, A_prod) * U)
            end
        end
    end

    # Rotate all observables in each unit into local reference frame.
    for unit in 1:nunits
        U = local_unitaries[unit]

        # Rotate interactions into local reference frames
        int = sys_contracted.interactions_union[unit]

        # Rotate onsite anisotropy (not that, for entangled units, onsite already includes Zeeman)
        int.onsite = Hermitian(U' * int.onsite * U)

        # Transform pair couplings into tensor decomposition and rotate.
        pair_new = PairCoupling[]
        for pc in int.pair
            # Convert PairCoupling to a purely general (tensor decomposed) interaction.
            pc_general = as_general_pair_coupling(pc, sys_contracted)

            # Rotate tensor decomposition into local frame.
            bond = pc.bond
            @assert bond.i == unit
            U′ = local_unitaries[bond.j]
            pc_rotated = rotate_general_coupling_into_local_frame(pc_general, U, U′)

            push!(pair_new, pc_rotated)
        end
        int.pair = pair_new
    end

    return SWTDataEntangled(
        local_unitaries,
        observables_localized,
        observable_buf,
    )
end


function intensities_bands(swt::EntangledSpinWaveTheory, qpts; kT=0)
    (; sys, contraction_info, data, measure, crystal_origin) = swt
    isempty(measure.observables) && error("No observables! Construct SpinWaveTheory with a `measure` argument.")

    qpts = convert(AbstractQPoints, qpts)
    rs_global = global_positions(sys)

    # Number of atoms in magnetic cell
    @assert sys.dims == (1,1,1)
    nunits = nsites(sys)                  # System "sites" are really entangled units
    nunits_per_cell = natoms(sys.crystal) # Crystal "atoms" are really entangled units
    # Number of chemical cells in magnetic cell
    Ncells = nunits / nunits_per_cell
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Preallocation
    T = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    Avec_pref = zeros(ComplexF64, nunits)
    disp = zeros(Float64, L, Nq)
    intensity = zeros(eltype(measure), L, Nq)

    # Temporary storage for pair correlations
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    Nobs = size(measure.observables, 1)

    for (iq, q) in enumerate(qpts.qs)
        q_global = crystal_origin.recipvecs * q
        q_reshaped = sys.crystal.recipvecs \ (crystal_origin.recipvecs * q)
        view(disp, :, iq) .= view(excitations!(T, H, swt, q), 1:L)

        for i in 1:nunits
            Avec_pref[i] = exp(- im * dot(q_global, rs_global[i]))
        end

        Avec = zeros(ComplexF64, Nobs)

        # Fill `intensity` array
        O = data.observable_buf
        for band = 1:L
            fill!(Avec, 0)
            N = sys.Ns[1]
            t = reshape(view(T, :, band), N-1, nunits, 2)
            for i in 1:nunits
                inverse_infos = contraction_info.inverse[i]
                for μ in 1:Nobs
                    fill!(O, 0)
                    for (subsite, inverse_info) in enumerate(inverse_infos)
                        site_original, offset = inverse_info.site, inverse_info.offset
                        ff = get_swt_formfactor(measure, 1, site_original)
                        prefactor = exp(-2π*im*(q_reshaped ⋅ offset)) * compute_form_factor(ff, norm2(q_global))
                        O += prefactor * data.observables_localized[subsite, μ, i]
                    end
                    for α in 1:N-1
                        Avec[μ] += Avec_pref[i] * (O[α, N] * t[α, i, 2] + O[N, α] * t[α, i, 1])
                    end
                end
            end

            map!(corrbuf, measure.corr_pairs) do (μ, ν)
                Avec[μ] * conj(Avec[ν]) / Ncells
            end
            intensity[band, iq] = thermal_prefactor(disp[band, iq]; kT) * measure.combiner(q_global, corrbuf)
        end
    end

    disp = reshape(disp, L, size(qpts.qs)...)
    intensity = reshape(intensity, L, size(qpts.qs)...)
    return BandIntensities(crystal_origin, qpts, disp, intensity)
end

################################################################################
# To be simplified, but requires editing existing code
################################################################################
# No changes
function energy_per_site_lswt_correction(swt::EntangledSpinWaveTheory; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    Natoms = natoms(sys.crystal)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # The uniform correction to the classical energy (trace of the (1,1)-block
    # of the spin-wave Hamiltonian)
    dynamical_matrix!(H, swt, zero(Vec3))
    δE₁ = -real(tr(view(H, 1:L, 1:L))) / 2Natoms

    # Integrate zero-point energy over the first Brillouin zone 𝐪 ∈ [0, 1]³ for
    # magnetic cell in reshaped RLU
    δE₂ = hcubature((0,0,0), (1,1,1); opts...) do q_reshaped
        dynamical_matrix!(H, swt, q_reshaped)
        ωs = bogoliubov!(V, H)
        return sum(view(ωs, 1:L)) / 2Natoms
    end

    # Error bars in δE₂[2] are discarded
    return δE₁ + δE₂[1]
end

# Only change is no Ewald
function swt_hamiltonian_SUN!(H::Matrix{ComplexF64}, swt::EntangledSpinWaveTheory, q_reshaped::Vec3)
    (; sys) = swt

    N = sys.Ns[1]
    Na = natoms(sys.crystal)
    L = (N-1) * Na

    # Clear the Hamiltonian
    @assert size(H) == (2L, 2L)
    H .= 0
    blockdims = (N-1, Na, N-1, Na)
    H11 = reshape(view(H, 1:L, 1:L), blockdims)
    H12 = reshape(view(H, 1:L, L+1:2L), blockdims)
    H21 = reshape(view(H, L+1:2L, 1:L), blockdims)
    H22 = reshape(view(H, L+1:2L, L+1:2L), blockdims)

    for (i, int) in enumerate(sys.interactions_union)

        # Onsite coupling, including Zeeman. Note that op has already been
        # transformed according to the local frame of sublattice i.
        op = int.onsite
        for m in 1:N-1
            for n in 1:N-1
                c = op[m, n] - δ(m, n) * op[N, N]
                H11[m, i, n, i] += c
                H22[n, i, m, i] += c
            end
        end

        for coupling in int.pair
            (; isculled, bond) = coupling
            isculled && break

            @assert i == bond.i
            j = bond.j

            phase = exp(2π*im * dot(q_reshaped, bond.n)) # Phase associated with periodic wrapping

            # Set "general" pair interactions of the form Aᵢ⊗Bⱼ. Note that Aᵢ
            # and Bᵢ have already been transformed according to the local frames
            # of sublattice i and j, respectively.
            for (Ai, Bj) in coupling.general.data
                for m in 1:N-1, n in 1:N-1
                    c = (Ai[m,n] - δ(m,n)*Ai[N,N]) * (Bj[N,N])
                    H11[m, i, n, i] += c
                    H22[n, i, m, i] += c

                    c = Ai[N,N] * (Bj[m,n] - δ(m,n)*Bj[N,N])
                    H11[m, j, n, j] += c
                    H22[n, j, m, j] += c

                    c = Ai[m,N] * Bj[N,n]
                    H11[m, i, n, j] += c * phase
                    H22[n, j, m, i] += c * conj(phase)

                    c = Ai[N,m] * Bj[n,N]
                    H11[n, j, m, i] += c * conj(phase)
                    H22[m, i, n, j] += c * phase

                    c = Ai[m,N] * Bj[n,N]
                    H12[m, i, n, j] += c * phase
                    H12[n, j, m, i] += c * conj(phase)
                    H21[n, j, m, i] += conj(c) * conj(phase)
                    H21[m, i, n, j] += conj(c) * phase
                end
            end
        end
    end

    # H must be hermitian up to round-off errors
    @assert maxdiff(H, H') < 1e-12

    # Make H exactly hermitian
    hermitianpart!(H)

    # Add small constant shift for positive-definiteness
    for i in 1:2L
        H[i,i] += swt.regularization
    end
end
