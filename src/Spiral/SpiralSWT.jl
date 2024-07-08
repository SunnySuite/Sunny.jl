
function construct_uniaxial_anisotropy(; axis, c20=0., c40=0., c60=0., S)
    # Anisotropy operator in local frame
    O = stevens_matrices(S)
    op = c20*O[2, 0] + c40*O[4, 0] + c60*O[6, 0]
    # Rotate operator into global frame, defined by axis
    R = Sunny.rotation_between_vectors(axis, [0, 0, 1])
    return rotate_operator(op, R)
end


## Dispersion and intensities

function swt_hamiltonian_dipole_spiral!(H::Matrix{ComplexF64}, swt::SpinWaveTheory, q_reshaped; k, axis)
    (; sys, data) = swt
    (; local_rotations, stevens_coefs) = data
    L = Sunny.nbands(swt) 
    @assert size(H) == (2L, 2L)
    H .= 0.0 
    
    # Add pairwise bilinear term
    for ints in sys.interactions_union

        for c in ints.pair
            (; i, j, n) = c.bond
            θ = (2*π * dot(k,n))
            Rn = axis_angle_to_matrix(axis, θ)

            # Undo rotations that were created in SpinWaveTheory.jl
            Ri = local_rotations[i]
            Rj = local_rotations[j]
            J = Ri * c.bilin * Rj'
            
            Jij = (J * Rn + Rn * J) ./ 2
            phase = exp(2π * im * dot(q_reshaped, n))
            
            Si = (sys.Ns[i]-1)/2
            Sj = (sys.Ns[j]-1)/2  

            ui = Ri[:,1] + im*Ri[:,2]
            uj = Rj[:,1] + im*Rj[:,2]
            vi = Ri[:,3]
            vj = Rj[:,3]
            
            H[i,j]     += (sqrt(Si*Sj)/2) * (transpose(ui)) * Jij * conj(uj) * phase
            H[i+L,j+L] += (sqrt(Si*Sj)/2) * conj((transpose(ui)) * Jij * conj(uj)) * phase
          
            H[i,j+L]   += (sqrt(Si*Sj)/2) * (transpose(ui) * Jij * uj) * phase
            H[j+L,i]   += (sqrt(Si*Sj)/2) * conj(transpose(ui) * Jij * uj * phase)
          
            H[i,i]     -= Sj * transpose(vi) * Jij * vj 
            H[i+L,i+L] -= Sj * transpose(vi) * Jij * vj

            iszero(c.biquad) || error("Biquadratic interactions not supported")
        end
    end

    H[:,:] = H / 2

    # Add Zeeman term
    for i in 1:L
        B = sys.extfield[1, 1, 1, i]' * sys.gs[1, 1, 1, i]
        B′ = - (B * local_rotations[i][:, 3]) / 2 
        H[i, i]     += B′
        H[i+L, i+L] += conj(B′)
    end
    
    # Add onsite couplings
    for i in 1:L
        S = (sys.Ns[i]-1)/2
        (; c2, c4, c6) = stevens_coefs[i]
        H[i, i]     += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i+L, i+L] += -3S*c2[3] - 40*S^3*c4[5] - 168*S^5*c6[7]
        H[i, i+L]   += +im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
        H[i+L, i]   += -im*(S*c2[5] + 6S^3*c4[7] + 16S^5*c6[9]) + (S*c2[1] + 6S^3*c4[3] + 16S^5*c6[5])
    end

    isnothing(sys.ewald) || error("Ewald interactions not yet supported")
        
    @assert diffnorm2(H, H') < 1e-12
    Sunny.hermitianpart!(H)
    
    for i in 1:2L
        H[i, i] += swt.energy_ϵ
    end
end

function dispersion_spiral(swt::SpinWaveTheory, axis; k, qs)
    (; sys) = swt
    
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm
   
    disp = zeros(Float64, nmodes, length(qs),3)

    for (iq, q) in enumerate(qs)
        for branch = 1:3    # 3 branch corresponds to K,K+Q and K-Q modes of incommensurate spin structures.
            H = zeros(ComplexF64, 2nmodes, 2nmodes)
            V = zeros(ComplexF64, 2nmodes, 2nmodes)
            q_reshaped = Sunny.to_reshaped_rlu(swt.sys, q)
            if sys.mode == :SUN
                error("Spiral calculation for SUN is not yet implemented")
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                swt_hamiltonian_dipole_spiral!(H, swt, q_reshaped .+ (branch - 2) .* k; k, axis)
            end
            try
                view(disp, :, iq,branch) .= Sunny.bogoliubov!(V, H)
            catch e
                error("Instability at wavevector q = $q")
            end
        end
    end

    return disp 
end


struct DipoleSpiralSpinWaveIntensityFormula{T}
    string_formula :: String
    kernel :: Union{Nothing,Function}
    calc_intensity :: Function
end

function Base.show(io::IO, ::DipoleSpiralSpinWaveIntensityFormula{T}) where T
    print(io,"SpinWaveIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::DipoleSpiralSpinWaveIntensityFormula{T}) where T
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
        println(io, "BandStructure information (ωᵢ and intensity) reported for each band")
    else
        println(io, "Intensity(ω) reported")
    end
end


function intensity_formula_spiral(swt::SpinWaveTheory, mode::Symbol; k, axis, kwargs...)
    contractor, string_formula = Sunny.contractor_from_mode(swt, mode)
    intensity_formula_spiral(swt, contractor; k, axis, string_formula, kwargs...)
end

function intensity_formula_spiral(swt::SpinWaveTheory, contractor::Sunny.Contraction{T}; k, axis, kwargs...) where T
    intensity_formula_spiral(swt, Sunny.required_correlations(contractor); k, axis, return_type=T, kwargs...) do qs, ωs, correlations
        intensity = Sunny.contract(correlations, qs, contractor)
    end
end

"""
    intensity_formula_spiral(swt, [contraction_mode]; k, axis)

Establish a formula for computing the scattering intensity by diagonalizing the
hamiltonian ``H(q)`` using Linear Spin Wave Theory.

The optional `contraction_mode` argument may be one of:
- `:trace` (default), which yields ``\\operatorname{tr} 𝒮(q,ω) = ∑_α
  𝒮^{αα}(q,ω)``
- `:perp`, which contracts ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} -
  q_{α}q_{β}``, returning the unpolarized intensity.
- `:full`, which will return all elements ``𝒮^{αβ}(𝐪,ω)`` without contraction.

If `kernel=delta_function_kernel`, then the resulting formula can be used with
[`intensities_bands`](@ref).

If `kernel` is an energy broadening kernel function, then the resulting formula
can be used with [`intensities_broadened`](@ref). Energy broadening kernel
functions can either be a function of `Δω` only, e.g.:

    kernel = Δω -> ...

or a function of both the energy transfer `ω` and of `Δω`, e.g.:

    kernel = (ω,Δω) -> ...

The integral of a properly normalized kernel function over all `Δω` is one.
"""
function intensity_formula_spiral(f::Function, swt::SpinWaveTheory, corr_ix::AbstractVector{Int64};
                           k, axis, kernel::Union{Nothing,Function},
                           return_type=Float64, string_formula="f(Q,ω,S{α,β}[ix_q,ix_ω])", 
                           formfactors=nothing)
    (; sys, data, observables) = swt
    Na = length(sys.dipoles) # number of magnetic atoms
    L = Sunny.nbands(swt) # k, k+Q, k-Q

    # Rotation matrices associated with `axis`
    CMat3 = SMatrix{3, 3, ComplexF64, 9}
    nx = CMat3([0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0])
    R2 = CMat3(axis * axis')
    R1 = (1/2) .* CMat3(I - im .* nx - R2)

    H = zeros(ComplexF64, 2L, 2L)
    T = zeros(ComplexF64, 2L, 2L, 3)
    tmp = zeros(ComplexF64, 2L, 2L)

    disp = zeros(Float64, L, 3)
    intensity = zeros(return_type, L,3)
    S = zeros(ComplexF64,3,3,L,3)

    FF = zeros(ComplexF64, Na)
    #intensity = zeros(return_type, nmodes,3)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = Sunny.propagate_form_factors_to_atoms(formfactors, swt.sys.crystal)
    
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
    calc_intensity = function(swt::SpinWaveTheory, q::Sunny.Vec3)
        # This function, calc_intensity, is an internal function to be stored
        # inside a formula. The unit system for `q` that is passed to
        # formula.calc_intensity is an implementation detail that may vary
        # according to the "type" of a formula. In the present context, namely
        # LSWT formulas, `q` is given in RLU for the original crystal. This
        # convention must be consistent with the usage in various
        # `intensities_*` functions defined in LinearSpinWaveIntensities.jl.
        # Separately, the functions calc_intensity for formulas associated with
        # SampledCorrelations will receive `q_absolute` in absolute units.
        
        q_reshaped = Sunny.to_reshaped_rlu(swt.sys, q)
        q_absolute = swt.sys.crystal.recipvecs * q_reshaped 

        for branch = 1:3   # 3 branch corresponds to K,K+Q and K-Q modes of incommensurate spin structures.
            if sys.mode == :SUN
                error("Spiral calculation for SUN is not yet implemented")
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                
                swt_hamiltonian_dipole_spiral!(H, swt, q_reshaped + (branch-2)*k; k, axis)
            
                disp[:,branch] = try
                    Sunny.bogoliubov!(tmp, H)
                catch e
                    error("Instability at wavevector q = $q")
                end

                T[:,:,branch] = tmp
            end
        end

        for i = 1:Na
            @assert Na == Sunny.natoms(sys.crystal)
            # TODO: move form factor into `f`, then delete this rescaling
            if isnothing(formfactors)
                FF[i] = 1.0
            else
                FF[i] = Sunny.compute_form_factor(ff_atoms[i], q_absolute⋅q_absolute)
            end
        end
       
        R = data.local_rotations
        Y = zeros(ComplexF64,L,L,3,3)
        Z = zeros(ComplexF64,L,L,3,3)
        V = zeros(ComplexF64,L,L,3,3)
        W = zeros(ComplexF64,L,L,3,3)
        for α in 1:3, β in 1:3
            for i in 1:L, j in 1:L
                si = (sys.Ns[i]-1)/2
                sj = (sys.Ns[j]-1)/2
                R_i = R[i]
                R_j = R[j]
                ui = R_i[:,1]+im*R_i[:,2]
                uj = R_j[:,1]+im*R_j[:,2]
                ti = sys.crystal.positions[i]
                tj = sys.crystal.positions[j]
                phase = exp(-2π * im*dot(q_reshaped, tj-ti))
                Y[i,j,α,β] = FF[i]*FF[j]*sqrt(si*sj) * (ui[α] * conj(uj[β])) * (phase)
                Z[i,j,α,β] = FF[i]*FF[j]*sqrt(si*sj) * (ui[α] * uj[β]) * (phase)
                V[i,j,α,β] = FF[i]*FF[j]*sqrt(si*sj) * (conj(ui[α]) * conj(uj[β])) * (phase)
                W[i,j,α,β] = FF[i]*FF[j]*sqrt(si*sj) * (conj(ui[α]) * uj[β]) * (phase)
            end
        end
        YZVW = [[Y Z];[V W]]

        for branch = 1:3, band = 1:L
            if sys.mode == :SUN
                error("Spiral calculation for SUN is not yet implemented")
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                for α in 1:3
                    for β in 1:3
                        A = T[:,:,branch]' * YZVW[:,:,α,β] * T[:,:,branch]
                        S[α,β,band,branch] = (1/(2*Na)) * A[band,band] 
                    end
                end
            end
        end
        
        avg(S) = 1/2 * (S - nx * S * nx + (R2 - I) * S * R2 + R2 * S * (R2 -I) + R2 * S * R2)
        
        for band = 1:L
            S[:,:,band,1] = avg(CMat3(S[:,:,band,1])) * conj(R1)
            S[:,:,band,2] = avg(CMat3(S[:,:,band,2])) * R2
            S[:,:,band,3] = avg(CMat3(S[:,:,band,3])) * R1
        end
        
        for branch = 1:3, band = 1:L
            @assert observables.observable_ixs[:Sx] == 1
            @assert observables.observable_ixs[:Sy] == 2
            @assert observables.observable_ixs[:Sz] == 3

            corrs = Vector{ComplexF64}(undef, Sunny.num_correlations(observables))
            for (ci,i) in observables.correlations
                (α,β) = ci.I

                corrs[i] = S[α,β,band,branch]
            end
            
            intensity[band, branch] = f(q_absolute, disp[band,branch], corrs[corr_ix])
        end
    
        # Return the result of the diagonalization in an appropriate
        # format based on the kernel provided
        if isnothing(kernel)
            # Delta function kernel --> (eigenvalue,intensity) pairs

            # If there is no specified kernel, we are done. Sort the bands in
            # order of decreasing dispersion, and return the BandStructure
            P = sortperm(vec(disp); rev=true)
            return Sunny.BandStructure{3*L,return_type}(disp[P], intensity[P])
        else
            disp_all = reshape(disp,:)
            intensity_all = reshape(intensity,:)
            # Smooth kernel --> Intensity as a function of ω (or a list of ωs)
            return function(ω)
                is = Vector{return_type}(undef,length(ω))
                is .= sum(intensity_all' .* kernel_edep.(disp_all', ω .- disp_all'),dims=2)
                is
            end
        end
    end
    output_type = isnothing(kernel) ? Sunny.BandStructure{L,return_type} : return_type
    DipoleSpiralSpinWaveIntensityFormula{output_type}(string_formula, kernel_edep, calc_intensity)
end

function intensities_bands(swt::SpinWaveTheory, qs, formula::DipoleSpiralSpinWaveIntensityFormula)
    if !isnothing(formula.kernel)
        # This is only triggered if the user has explicitly specified a formula with e.g. kT
        # corrections applied, but has not disabled the broadening kernel.
        error("intensities_bands: Can't compute band intensities if a broadening kernel is applied.\nTry intensity_formula(...; kernel = delta_function_kernel)")
    end

    qs = Sunny.Vec3.(qs)
    nmodes = Sunny.nbands(swt)

    # Get the type parameter from the BandStructure
    return_type = typeof(formula).parameters[1].parameters[2]

    band_dispersions = zeros(Float64, size(qs)..., 3*nmodes)
    band_intensities = zeros(return_type, size(qs)..., 3*nmodes)
    for qidx in CartesianIndices(qs)
        band_structure = formula.calc_intensity(swt, qs[qidx])

        band_dispersions[qidx,:] .= band_structure.dispersion
        band_intensities[qidx,:] .= band_structure.intensity
    end
    return band_dispersions, band_intensities
end
