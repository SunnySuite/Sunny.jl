###########################################################################
# Below takes Sunny to construct `SpinWave` for LSWT calculations.  #
###########################################################################
struct SWTDataDipole
    R_mat  :: Vector{Mat3}             # SO(3) rotation to align the quantization axis
    c_coef :: Vector{StevensExpansion} # Stevens operator coefficents
end

struct SWTDataSUN
    dipole_operators        :: Array{ComplexF64, 4}
    quadrupole_operators    :: Array{ComplexF64, 4}
    onsite_operator         :: Array{ComplexF64, 3}  
    external_field_operator :: Array{ComplexF64, 3}
    general_pair_operators  :: Vector{Tuple{Tuple{HermitianC64, HermitianC64}, Bond}}
    observable_operators    :: Array{ComplexF64, 4}
end

"""
    SpinWaveTheory(sys, energy_ϵ::Float64=1e-8, energy_tol=1e-6)

Constructs an object to perform linear spin wave theory. Use it with
[`dispersion`](@ref) and [`dssf`](@ref) functions.

The optional parameter `energy_ϵ` adds a small positive shift to the diagonal of
the dynamical matrix ``D`` to avoid numerical issues with zero-energy
quasi-particle modes. The optional parameter `energy_tol` relaxes the check on
the imaginary part of the eigenvalues.
"""
struct SpinWaveTheory
    sys        :: System
    data       :: Union{SWTDataDipole, SWTDataSUN}
    energy_ϵ   :: Float64          # Energy shift applied to dynamical matrix prior to Bogoliubov transformation
    energy_tol :: Float64          # Energy tolerance for maximal imaginary part of quasiparticle energies

    observables :: ObservableInfo
end

function SpinWaveTheory(sys::System{N}; energy_ϵ::Float64=1e-8, energy_tol::Float64=1e-6, observables=nothing, correlations=nothing) where N
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
        data = swt_data_dipole(sys)
    end

    return SpinWaveTheory(sys, data, energy_ϵ, energy_tol, obs)
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


# Convert 3-vector from the Cartesian frame to the spherical frame
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

# Compute SU(N) generators in the local reference frame (for :SUN mode).
function swt_data_sun(sys::System{N}, obs) where N
    n_magnetic_atoms = natoms(sys.crystal)

    dipole_operators = spin_matrices_of_dim(; N)
    quadrupole_operators = stevens_matrices_of_dim(2; N)

    # Rotate dipole and quadrupoles into local reference frames
    local_quantization_bases = [Matrix{ComplexF64}(undef, N, N) for _ in 1:n_magnetic_atoms]

    for atom in 1:n_magnetic_atoms
        # First axis of local quantization basis is along the 
        # ground-state polarization axis
        local_quantization_bases[atom][:, 1] = sys.coherents[1,1,1,atom]
        
        # Remaining axes are arbitrary but mutually orthogonal
        # and orthogonal to the first axis
        local_quantization_bases[atom][:, 2:N] = nullspace(local_quantization_bases[atom][:,1]') 
    end

    dipole_operators_localized = zeros(ComplexF64, 3, N, N, n_magnetic_atoms)
    onsite_operator_localized = zeros(ComplexF64, N, N, n_magnetic_atoms)
    quadrupole_operators_localized = zeros(ComplexF64, 5, N, N, n_magnetic_atoms)
    observables_localized = zeros(ComplexF64, N, N, num_observables(obs), n_magnetic_atoms)

    # Rotate SU(N) bases and observables and store in dense array. Note that the
    # first index is the component index. As a result, (Sˣᵢⱼ, Sʸᵢⱼ, Sᶻᵢⱼ) is
    # stored contiguously for each matrix element ij. This is the natural order
    # for constructing the spin wave Hamiltonian.
    for atom in 1:n_magnetic_atoms
        U = local_quantization_bases[atom]
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
    # frames.
    (; extfield, gs, units) = sys
    external_field_operator = zeros(ComplexF64, N, N, n_magnetic_atoms) 
    for atom in 1:n_magnetic_atoms
        B = units.μB * (gs[1, 1, 1, atom]' * extfield[1, 1, 1, atom])
        S = view(dipole_operators_localized, :, :, :, atom)
        @. external_field_operator[:, :, atom] = -B[1]*S[1, :, :] - B[2]*S[2, :, :] - B[3]*S[3, :, :]
    end

    # Rotate operators generated by the tensor decomposition of generalized
    # interactions.
    general_pair_operators_localized = Tuple{Tuple{HermitianC64, HermitianC64}, Bond}[]
    for int in sys.interactions_union
        for pc in int.pair
            (; isculled, bond, general) = pc
            isculled && continue 
            atomi, atomj = bond.i, bond.j
            Ui, Uj = local_quantization_bases[atomi], local_quantization_bases[atomj] 
            for (A, B) in general.data
                A′, B′ = Hermitian.((Ui' * A * Ui, Uj' * B * Uj)) 
                push!(general_pair_operators_localized, ((A′, B′), bond))
            end
        end
    end

    return SWTDataSUN(
        dipole_operators_localized, 
        quadrupole_operators_localized, 
        onsite_operator_localized,
        external_field_operator,
        general_pair_operators_localized,
        observables_localized,
    )
end

# Compute Stevens coefficients in the local reference frame
function swt_data_dipole(sys::System{0})
    cs = StevensExpansion[]
    Rs = Mat3[]

    for atom in 1:natoms(sys.crystal)
        # SO(3) rotation that aligns the quantization axis. Important: since we
        # will project out bosons that correspond to multipolar fluctuations,
        # therefore we use the explicit matrix to get rid of any ambiguity.
        #
        # As a unitary, U = exp(-i ϕ Sz) exp(-i θ Sy)
        θ, ϕ = dipole_to_angles(sys.dipoles[1,1,1,atom])
        R = SA[-sin(ϕ) -cos(ϕ)*cos(θ) cos(ϕ)*sin(θ);
                cos(ϕ) -sin(ϕ)*cos(θ) sin(ϕ)*sin(θ);
                0.0     sin(θ)        cos(θ)]
        # Rotated Stevens expansion
        c = rotate_operator(sys.interactions_union[atom].onsite, R)

        push!(Rs, R)
        push!(cs, c)
    end

    return SWTDataDipole(Rs, cs)
end
