# Bogoliubov transformation that diagonalizes a quadratic bosonic Hamiltonian,
# allowing for anomalous terms. The general procedure derives from Colpa,
# Physica A, 93A, 327-353 (1978).
function bogoliubov!(V::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    L = div(size(H, 1), 2)
    @assert size(V) == size(H) == (2L, 2L)

    # Initialize V to the para-unitary identity Ĩ = diagm([ones(L), -ones(L)])
    V .= 0
    for i in 1:L
        V[i, i] = 1
        V[i+L, i+L] = -1
    end

    # Solve generalized eigenvalue problem, Ĩ t = λ H t, for columns t of V.
    # Eigenvalues are sorted such that quasi-particle energies will appear in
    # descending order.
    λ, V0 = eigen!(Hermitian(V), Hermitian(H); sortby = x -> -1/real(x))

    # Note that V0 and V refer to the same data.
    @assert V0 === V

    # Normalize columns of V so that para-unitarity holds, V† Ĩ V = Ĩ.
    for j in axes(V, 2)
        c = 1 / sqrt(abs(λ[j]))
        view(V, :, j) .*= c
    end

    # Disable test for speed
    #=
    Ĩ = Diagonal([ones(L); -ones(L)])
    @assert V' * Ĩ * V ≈ Ĩ
    =#

    # Verify that half the eigenvalues are positive and the other half are
    # negative. The positive eigenvalues are quasiparticle energies for the
    # wavevector q that defines the dynamical matrix H(q). The absolute value of
    # the negative eigenvalues would be quasiparticle energies for H(-q), which
    # we are not considering in the present context.
    @assert all(>(0), view(λ, 1:L)) && all(<(0), view(λ, L+1:2L))
    
    # Inverse of λ gives eigenvalues of Ĩ H. We only care about the first L
    # eigenvalues, which are positive. A factor of 2 is needed to get the
    # physical quasiparticle energies.
    disp = resize!(λ, L)
    @. disp = 2 / disp

    # In the special case that H(q) = H(-q) (i.e., a magnetic ordering with
    # reflection symmetry), the eigenvalues come in pairs. Note that the data in
    # H has been overwritten by eigen!, so H0 should refer to an original copy
    # of H.
    #=
    @assert diag(V' * H0 * V) ≈ [disp/2; reverse(disp)/2]
    =#

    return disp
end


# DD: These two functions are a stopgap until data is treated differently in
# main calculations. Also, the final data layout will need to be iterated on. I
# am thinking the user should always be able to get some array with indices
# identical to the list of wave vectors. This could be achieved, for example, by
# having the output be an array with length equal to the number of modes. Each
# entry would then be an array with dimension equal to the array of wave
# vectors. The entries of this array would then depend on the request (an an
# energy, an intensity, an entire tensor stored as an SMatrix, etc.) 
# The key point is to make it as easy as possible to put the output
# in correspondence with the input for plotting, further processing, etc.
function reshape_correlations(corrs)
    qdims, nmodes = size(corrs)[4:end], size(corrs)[3]  # First two indices are are always tensor indices
    idxorder = collect(1:ndims(corrs))
    idxorder[3], idxorder[end] = idxorder[end], idxorder[3]
    corrs = permutedims(corrs, idxorder)
    return selectdim(reinterpret(SMatrix{3,3,ComplexF64,9}, reshape(corrs, 9, qdims...,nmodes) ), 1, 1)
end

function reshape_dispersions(disp)
    idxorder = collect(1:ndims(disp))
    idxorder[1], idxorder[end] = idxorder[end], idxorder[1]
    return permutedims(disp, idxorder)
end

"""
    dispersion(swt::SpinWaveTheory, qs)

Computes the spin excitation energy dispersion relations given a
[`SpinWaveTheory`](@ref) and an array of wave vectors `qs`. Each element ``q``
of `qs` must be a 3-vector in units of reciprocal lattice units. I.e., ``qᵢ`` is
given in ``2π/|aᵢ|`` with ``|aᵢ|`` the lattice constant of the original chemical
lattice.

The first indices of the returned array correspond to those of `qs`. A final
index, corresponding to mode, is added to these. Each entry of the array is an
energy.
"""
function dispersion(swt::SpinWaveTheory, qs)
    (; sys) = swt
    
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm

    H = zeros(ComplexF64, 2nmodes, 2nmodes)
    V = zeros(ComplexF64, 2nmodes, 2nmodes)
    disp = zeros(Float64, nmodes, length(qs))

    for (iq, q) in enumerate(qs)
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        if sys.mode == :SUN
            swt_hamiltonian_SUN!(H, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H, swt, q_reshaped)
        end

        try
            view(disp, :, iq) .= bogoliubov!(V, H)
        catch e
            error("Instability at wavevector q = $q")
        end
    end

    return reshape_dispersions(disp)
end


"""
    dssf(swt::SpinWaveTheory, qs)

Given a [`SpinWaveTheory`](@ref) object, computes the dynamical spin structure
factor,

```math
    𝒮^{αβ}(𝐤, ω) = 1/(2πN)∫dt ∑_𝐫 \\exp[i(ωt - 𝐤⋅𝐫)] ⟨S^α(𝐫, t)S^β(0, 0)⟩,
```

using the result from linear spin-wave theory,

```math
    𝒮^{αβ}(𝐤, ω) = ∑_n |A_n^{αβ}(𝐤)|^2 δ[ω-ω_n(𝐤)].
```

`qs` is an array of wave vectors of arbitrary dimension. Each element ``q`` of
`qs` must be a 3-vector in reciprocal lattice units (RLU), i.e., in the basis of
reciprocal lattice vectors.

The first indices of the returned array correspond to those of `qs`. A final
index, corresponding to mode, is added to these. Each entry of this array is a
tensor (3×3 matrix) corresponding to the indices ``α`` and ``β``.
"""
function dssf(swt::SpinWaveTheory, qs)
    qs = Vec3.(qs)
    nmodes = nbands(swt)

    disp = zeros(Float64, nmodes, size(qs)...)
    Sαβs = zeros(ComplexF64, 3, 3, nmodes, size(qs)...) 

    # dssf(...) doesn't do any contraction, temperature correction, etc.
    # It simply returns the full Sαβ correlation matrix
    formula = intensity_formula(swt, :full; kernel = delta_function_kernel)

    # Calculate DSSF 
    for qidx in CartesianIndices(qs)
        q = qs[qidx]
        band_structure = formula.calc_intensity(swt,q)
        for band = 1:nmodes
            disp[band,qidx] = band_structure.dispersion[band]
            Sαβs[:,:,band,qidx] .= band_structure.intensity[band]
        end
    end

    return reshape_dispersions(disp), reshape_correlations(Sαβs) 
end 


struct BandStructure{N,T}
  dispersion :: SVector{N,Float64}
  intensity :: SVector{N,T}
end

struct SpinWaveIntensityFormula{T}
    string_formula :: String
    kernel :: Union{Nothing,Function}
    calc_intensity :: Function
end

function Base.show(io::IO, ::SpinWaveIntensityFormula{T}) where T
    print(io,"SpinWaveIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::SpinWaveIntensityFormula{T}) where T
    printstyled(io, "Quantum Scattering Intensity Formula\n"; bold=true, color=:underline)

    formula_lines = split(formula.string_formula, '\n')

    if isnothing(formula.kernel)
        println(io, "At any Q and for each band ωᵢ = εᵢ(Q), with S = S(Q,ωᵢ):\n")
        intensity_equals = "  Intensity(Q,ω) = ∑ᵢ δ(ω-ωᵢ) "
    else
        println(io, "At any (Q,ω), with S = S(Q,ωᵢ):\n")
        intensity_equals = "  Intensity(Q,ω) = ∑ᵢ Kernel(ω-ωᵢ) "
    end
    separator = '\n' * repeat(' ', textwidth(intensity_equals))
    println(io, intensity_equals, join(formula_lines, separator))
    println(io)
    if isnothing(formula.kernel)
        println(io,"BandStructure information (ωᵢ and intensity) reported for each band")
    else
        println(io,"Intensity(ω) reported")
    end
end

delta_function_kernel = nothing

"""
    formula = intensity_formula(swt::SpinWaveTheory; kernel = ...)

Establish a formula for computing the scattering intensity by diagonalizing
the hamiltonian ``H(q)`` using Linear Spin Wave Theory.

If `kernel = delta_function_kernel`, then the resulting formula can be used with
[`intensities_bands`](@ref).

If `kernel` is an energy broadening kernel function, then the resulting formula can be used with [`intensities_broadened`](@ref).
Energy broadening kernel functions can either be a function of `Δω` only, e.g.:

    kernel = Δω -> ...

or a function of both the energy transfer `ω` and of `Δω`, e.g.:

    kernel = (ω,Δω) -> ...

The integral of a properly normalized kernel function over all `Δω` is one.
"""
function intensity_formula(f::Function,swt::SpinWaveTheory,corr_ix::AbstractVector{Int64}; kernel::Union{Nothing,Function},
                           return_type=Float64, string_formula="f(Q,ω,S{α,β}[ix_q,ix_ω])", mode_fast=false,
                           formfactors=nothing)
    (; sys, data, observables) = swt
    Nm, N = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    S = (N-1) / 2
    nmodes = nbands(swt)
    sqrt_Nm_inv = 1.0 / √Nm
    sqrt_halfS  = √(S/2)

    # Preallocation
    H = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    V = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Avec_pref = zeros(ComplexF64, Nm)
    intensity = zeros(return_type, nmodes)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, swt.sys.crystal)

    # Upgrade to 2-argument kernel if needed
    kernel_edep = if isnothing(kernel)
        nothing
    else
        try
            kernel(0.,0.)
            kernel
        catch MethodError
            (ω,Δω) -> kernel(Δω)
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
    calc_intensity = function(swt::SpinWaveTheory, q::Vec3)
        # This function, calc_intensity, is an internal function to be stored
        # inside a formula. The unit system for `q` that is passed to
        # formula.calc_intensity is an implementation detail that may vary
        # according to the "type" of a formula. In the present context, namely
        # LSWT formulas, `q` is given in RLU for the original crystal. This
        # convention must be consistent with the usage in various
        # `intensities_*` functions defined in LinearSpinWaveIntensities.jl.
        # Separately, the functions calc_intensity for formulas associated with
        # SampledCorrelations will receive `q_absolute` in absolute units.
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        q_absolute = swt.sys.crystal.recipvecs * q_reshaped

        if sys.mode == :SUN
            swt_hamiltonian_SUN!(H, swt, q_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H, swt, q_reshaped)
        end

        disp = try
            bogoliubov!(V, H)
        catch e
            error("Instability at wavevector q = $q")
        end

        for i = 1:Nm
            @assert Nm == natoms(sys.crystal)
            phase = exp(-2π*im * dot(q_reshaped, sys.crystal.positions[i]))
            Avec_pref[i] = sqrt_Nm_inv * phase

            # TODO: move form factor into `f`, then delete this rescaling
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)
        end

        # Fill `intensity` array
        for band = 1:nmodes
            corrs = if sys.mode == :SUN
                v = reshape(view(V, :, band), N-1, Nm, 2)
                Avec = zeros(ComplexF64, num_observables(observables))
                (; observable_operators) = data
                for i = 1:Nm
                    for μ = 1:num_observables(observables)
                        @views O = observable_operators[:, :, μ, i]
                        for α = 1:N-1
                            Avec[μ] += Avec_pref[i] * (O[α, N] * v[α, i, 2] + O[N, α] * v[α, i, 1])
                        end
                    end
                end
                corrs = Vector{ComplexF64}(undef,num_correlations(observables))
                for (ci,i) in observables.correlations
                    (α,β) = ci.I
                    corrs[i] = Avec[α] * conj(Avec[β])
                end
                corrs
            else
                v = reshape(view(V, :, band), Nm, 2)
                @assert sys.mode in (:dipole, :dipole_large_S)
                Avec = zeros(ComplexF64, 3)
                R = data.local_rotations
                for i = 1:Nm
                    Vtmp = [v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                    Avec += Avec_pref[i] * sqrt_halfS * (R[i] * Vtmp)
                end

                @assert observables.observable_ixs[:Sx] == 1
                @assert observables.observable_ixs[:Sy] == 2
                @assert observables.observable_ixs[:Sz] == 3
                corrs = Vector{ComplexF64}(undef, num_correlations(observables))
                for (ci, i) in observables.correlations
                    (α, β) = ci.I
                    corrs[i] = Avec[α] * conj(Avec[β])
                end
                corrs
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
                is .= sum(intensity' .* kernel_edep.(disp', ω .- disp'),dims=2)
                is
            end
        end
    end
    output_type = isnothing(kernel) ? BandStructure{nmodes,return_type} : return_type
    SpinWaveIntensityFormula{output_type}(string_formula, kernel_edep, calc_intensity)
end
