struct SWTDataEntangled
    local_unitaries           :: Array{ComplexF64, 3}  # Aligns to quantization axis on each site
    observable_operators_full :: Array{ComplexF64, 5}  # Observables in local frame for each subsite (for intensity calcs)
    observable_buf            :: Array{ComplexF64, 2}  # Buffer for use while constructing boson rep of observables 
end

struct EntangledSpinWaveTheory # Could just expand union above, but type now available for dispatch
    sys              :: System
    crystal_origin   :: Crystal
    contraction_info :: CrystalContractionInfo
    Ns_unit          :: Vector{Vector{Int64}}
    data             :: SWTDataEntangled
    energy_ϵ         :: Float64
    observables      :: ObservableInfo
end

function EntangledSpinWaveTheory(esys::EntangledSystem; energy_ϵ::Float64=1e-8, observables=nothing, correlations=nothing, apply_g = true)
    (; sys, sys_origin) = esys
    crystal_origin = orig_crystal(sys_origin)
    if !isnothing(sys.ewald)
        error("SpinWaveTheory does not yet support long-range dipole-dipole interactions.")
    end

    # Could move this to parse_observables, but then parse_observables would have to know about EntanglementData
    if isnothing(observables)
        # Note that the observables must have dimensions consistent with the original system
        Ns_original = vcat(Ns_in_units(sys_origin, esys.contraction_info)...)
        @assert allequal(Ns_original) "All local Hilbert spaces of original system must have identical dimension"
        N_original = Ns_original[1]
        S = spin_matrices((N_original-1)/2)
        observables = [
            :Sx => S[1],
            :Sy => S[2],
            :Sz => S[3],
        ]
    end

    # Reshape the entangled system so that it is a single unit cell.
    cellsize_mag = cell_shape(sys) * diagm(Vec3(sys.latsize))
    sys_reshaped = reshape_supercell_aux(sys, (1,1,1), cellsize_mag)

    # Reshape the original system in a corresponding fashion
    cellsize_original = cell_shape(sys_origin) * diagm(Vec3(sys_origin.latsize))
    sys_origin_reshaped = reshape_supercell_aux(sys_origin, (1,1,1), cellsize_original)

    # Construct the new contraction_info (i.e. reconstruct maps between
    # entangled and unentangled systems)
    new_units = units_for_reshaped_system(sys_origin_reshaped, esys)
    _, new_contraction_info = contract_crystal(sys_origin_reshaped.crystal, new_units)
    new_Ns_unit = Ns_in_units(sys_origin_reshaped, new_contraction_info)

    # Rotate local operators to quantization axis
    N = esys.sys.Ns[1]
    obs = parse_observables(N; observables, correlations, g = apply_g ? sys.gs : nothing)
    data = swt_data_entangled(sys_reshaped, new_Ns_unit, new_contraction_info, obs)

    return EntangledSpinWaveTheory(sys_reshaped, crystal_origin, new_contraction_info, new_Ns_unit, data, energy_ϵ, obs)
end

function Base.show(io::IO, ::MIME"text/plain", swt::EntangledSpinWaveTheory)
    printstyled(io, "SpinWaveTheory\n"; bold=true, color=:underline)
    println(io, "Entangled units in magnetic supercell: $(natoms(swt.sys.crystal))")
    show(io, MIME("text/plain"), swt.observables)
end

nbands(swt::EntangledSpinWaveTheory) = (swt.sys.Ns[1]-1)  * natoms(swt.sys.crystal)

function to_reshaped_rlu(esys::EntangledSystem, q)
    (; sys, sys_origin) = esys
    return sys.crystal.recipvecs \ (orig_crystal(sys_origin).recipvecs * q)
end

# obs are observables _given in terms of `sys_original`_
function swt_data_entangled(sys::System, Ns_unit, contraction_info, obs)
    N = sys.Ns[1]
    
    # Calculate transformation matrices into local reference frames
    nunits = natoms(sys.crystal)

    # Check to make sure all "units" are the same -- this can be generalized later.
    sites_per_unit = [length(info) for info in contraction_info.inverse]
    @assert allequal(sites_per_unit) "All units must have the same number of interior sites"
    sites_per_unit = sites_per_unit[1]

    Ns_contracted = map(Ns -> prod(Ns), Ns_unit)
    @assert allequal(Ns_contracted) "All units must have the same dimension local Hilbert space"
    @assert Ns_contracted[1] == N "Unit dimension inconsistent with system"  # Sanity check. This should never happen. 

    # Preallocate buffers for local unitaries and observables.
    local_unitaries = zeros(ComplexF64, N, N, nunits)
    observables_localized_all = zeros(ComplexF64, N, N, sites_per_unit, num_observables(obs), nunits)
    observable_buf = zeros(ComplexF64, N, N)

    for atom in 1:nunits
        # Create unitary that rotates [0, ..., 0, 1] into ground state direction
        # Z that defines quantization axis
        Z = sys.coherents[atom]
        view(local_unitaries, :, N, atom)     .= Z
        view(local_unitaries, :, 1:N-1, atom) .= nullspace(Z')
    end

    for unit in 1:nunits
        U = view(local_unitaries, :, :, unit)

        # Rotate observables into local reference frames
        for k = 1:num_observables(obs)
            A = obs.observables[k]
            for local_site in 1:sites_per_unit
                A_prod = local_op_to_product_space(A, local_site, Ns_unit[unit])
                observables_localized_all[:, :, local_site, k, unit] = Hermitian(U' * convert(Matrix, A_prod) * U)
            end
        end

        # Rotate interactions into local reference frames
        int = sys.interactions_union[unit]

        # Rotate onsite anisotropy (not that, for entangled units, onsite already includes Zeeman)
        int.onsite = Hermitian(U' * int.onsite * U) 

        # Transform pair couplings into tensor decomposition and rotate.
        pair_new = PairCoupling[]
        for pc in int.pair
            # Convert PairCoupling to a purely general (tensor decomposed) interaction.
            pc_general = as_general_pair_coupling(pc, sys)

            # Rotate tensor decomposition into local frame.
            bond = pc.bond
            @assert bond.i == unit
            U′ = view(local_unitaries, :, :, bond.j)
            pc_rotated = rotate_general_coupling_into_local_frame(pc_general, U, U′)

            push!(pair_new, pc_rotated)
        end
        int.pair = pair_new
    end

    return SWTDataEntangled(
        local_unitaries,
        observables_localized_all,
        observable_buf,
    )
end

function intensity_formula(swt::EntangledSpinWaveTheory, mode::Symbol; kwargs...)
    contractor, string_formula = contractor_from_mode(swt, mode)
    intensity_formula(swt, contractor; string_formula, kwargs...)
end

function intensity_formula(swt::EntangledSpinWaveTheory, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(swt, required_correlations(contractor); return_type=T, kwargs...) do k, _, correlations
        intensity = contract(correlations, k, contractor)
    end
end

function intensity_formula(f::Function, swt::EntangledSpinWaveTheory, corr_ix::AbstractVector{Int64}; 
    kernel::Union{Nothing, Function},
    return_type=Float64, 
    string_formula="f(Q,ω,S{α,β}[ix_q,ix_ω])", 
    formfactors=nothing  # TODO: Add support
)
    (; sys, contraction_info, data, observables) = swt

    nunits, N = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    nmodes = nunits * (N-1) 
    sqrt_nunits_inv = 1.0 / √nunits

    # Preallocation
    H = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    V = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Avec_pref = zeros(ComplexF64, nunits)
    intensity = zeros(return_type, nmodes)

    # TODO: form factors for entangled unit systems
    ff_atoms = propagate_form_factors_to_atoms(formfactors, swt.crystal_origin)

    # Upgrade to 2-argument kernel if needed
    kernel_edep = if isnothing(kernel)
        nothing
    else
        try
            kernel(0., 0.)
            kernel
        catch MethodError
            (ω, Δω) -> kernel(Δω)
        end
    end

    # In Spin Wave Theory, the Hamiltonian depends on momentum transfer `q`.
    # At each `q`, the Hamiltonian is diagonalized one time, and then the
    # energy eigenvalues can be reused multiple times. To facilitate this,
    # `I_of_ω = calc_intensity(swt,q)` performs the diagonalization, and returns
    # the result either as:
    #
    #   Delta function kernel --> I_of_ω = (eigenvalue,intensity) pairs
    #
    #   OR
    #
    #   Smooth kernel --> I_of_ω = Intensity as a function of ω
    #
    calc_intensity = function(swt::EntangledSpinWaveTheory, q::Vec3)
        # This function, calc_intensity, is an internal function to be stored
        # inside a formula. The unit system for `q` that is passed to
        # formula.calc_intensity is an implementation detail that may vary
        # according to the "type" of a formula. In the present context, namely
        # LSWT formulas, `q` is given in RLU for the original crystal. This
        # convention must be consistent with the usage in various
        # `intensities_*` functions defined in LinearSpinWaveIntensities.jl.
        # Separately, the functions calc_intensity for formulas associated with
        # SampledCorrelations will receive `q_absolute` in absolute units.
        (; sys, contraction_info, crystal_origin) = swt
        q_reshaped = sys.crystal.recipvecs \ (crystal_origin.recipvecs * q)
        q_absolute = sys.crystal.recipvecs * q_reshaped

        swt_hamiltonian_SUN!(H, swt, q_reshaped)

        disp = try
            bogoliubov!(V, H)
        catch e
            error("Instability at wavevector q = $q") end

        for i = 1:nunits
            @assert nunits == natoms(sys.crystal)
            phase = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))

            Avec_pref[i] = sqrt_nunits_inv * phase

            # TODO: FormFactors for entangled units 
            # Must be moved to where there is information about subsite
            # Avec_pref[i] *= compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)  
        end

        # Fill `intensity` array
        for band = 1:nmodes
            v = reshape(view(V, :, band), N-1, nunits, 2)
            Avec = zeros(ComplexF64, num_observables(observables))
            (; observable_buf, observable_operators_full) = data
            for i = 1:nunits
                inverse_infos = contraction_info.inverse[i]
                for μ = 1:num_observables(observables)

                    # Construct q-dependent observable
                    O = observable_buf
                    O .= 0
                    for (subsite, inverse_info) in enumerate(inverse_infos)
                        site_original, offset = inverse_info.site, inverse_info.offset
                        prefactor = exp(-2π*im*(q_reshaped ⋅ offset)) * compute_form_factor(ff_atoms[site_original], q_absolute ⋅ q_absolute)
                        O += prefactor * @view(observable_operators_full[:, :, subsite, μ, i])
                    end
                    unit_normalization = 1/√(length(inverse_infos))

                    # Accumulate into bosonic representation of observable
                    for α = 1:N-1
                        Avec[μ] += Avec_pref[i] * unit_normalization * (O[α, N] * v[α, i, 2] + O[N, α] * v[α, i, 1])
                    end
                end
            end

            # Calculate correlations
            corrs = Vector{ComplexF64}(undef,num_correlations(observables))
            for (ci, i) in observables.correlations
                (α,β) = ci.I
                corrs[i] = Avec[α] * conj(Avec[β])
            end

            intensity[band] = f(q_absolute, disp[band], corrs[corr_ix])
        end

        # Return the result of the diagonalization in an appropriate
        # format based on the kernel provided
        if isnothing(kernel)
            # Delta function kernel --> (eigenvalue,intensity) pairs

            # If there is no specified kernel, we are done: just return the
            # BandStructure
            return BandStructure{nmodes,return_type}(disp, intensity)
        else
            # Smooth kernel --> Intensity as a function of ω (or a list of ωs)
            return function(ω)
                is = Vector{return_type}(undef,length(ω))
                is .= sum(intensity' .* kernel_edep.(disp', ω .- disp'), dims=2)
                is
            end
        end
    end
    output_type = isnothing(kernel) ? BandStructure{nmodes,return_type} : return_type
    SpinWaveIntensityFormula{output_type}(string_formula, kernel_edep, calc_intensity)
end