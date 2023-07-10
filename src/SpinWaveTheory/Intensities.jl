#DD this will become more similar to the existing intensities.
"""
    intensities(swt::SpinWaveTheory, qs, ωvals, η::Float64)

Computes the unpolarized inelastic neutron scattering intensities given a
`SpinWaveTheory`, an array of wave vectors `qs`, a list of energies `ωvals`, and
a Lorentzian broadening parameter `η`.

Note that `qs` is an array of wave vectors of arbitrary dimension. Each element
``q`` of `qs` must be a 3-vector in reciprocal lattice units. I.e., ``qᵢ`` is
given in ``2π/|aᵢ|`` with ``|aᵢ|`` the lattice constant of the chemical lattice.

The output will be an array with indices identical to `qs`. Each entry of the
array will be an unpolarized intensity.
"""
# DD: incorporate existing SF utilties (e.g., form factor)
function intensities(swt::SpinWaveTheory, qs, ωvals, η::Float64; formula = intensity_formula(swt) :: SpinWaveIntensityFormula)
    (; sys) = swt
    qs = Vec3.(qs)
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    nmodes  = Nf * Nm

    num_ω = length(ωvals)
    is = zeros(Float64, size(qs)..., num_ω)
    if typeof(formula).parameters[1] != Float64
      error("Intensity formulae with complicated return types are not yet supported for SpinWaveTheory")
    end

    for qidx in CartesianIndices(qs)
        disp, intensity = formula.calc_intensity(swt,qs[qidx])

        for band = 1:nmodes
            # At a Goldstone mode, where the intensity is divergent, use a delta-function for the intensity.
            if (disp[band] < 1.0e-3) && (intensity[band] > 1.0e3)
                is[qidx, 1] += intensity[band]
            else
                for index_ω = 1:num_ω
                    is[qidx, index_ω] += intensity[band] * lorentzian(ωvals[index_ω]-disp[band], η)
                end
            end
        end
    end
    return is
end

struct SpinWaveIntensityFormula{T}
    kT :: Float64
    formfactors
    string_formula :: String
    calc_intensity :: Function
end

function Base.show(io::IO, formula::SpinWaveIntensityFormula{T}) where T
    print(io,"SpinWaveIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::SpinWaveIntensityFormula{T}) where T
    printstyled(io, "Quantum Scattering Intensity Formula\n";bold=true, color=:underline)

    formula_lines = split(formula.string_formula,'\n')

    intensity_equals = "  Intensity(Q,ω) = ∑ᵢ δ(ω-ωᵢ) "
    println(io,"At any Q and for each band ωᵢ = εᵢ(Q), with S = S(Q,ωᵢ):")
    println(io)
    println(io,intensity_equals,formula_lines[1])
    for i = 2:length(formula_lines)
        precursor = repeat(' ', textwidth(intensity_equals))
        println(io,precursor,formula_lines[i])
    end
    println(io)
    println(io,"Intensity :: $(T) and dispersion reported for each band")
end

intensity_formula(swt::SpinWaveTheory; kwargs...) = intensity_formula(swt, :perp; kwargs...)
function intensity_formula(swt::SpinWaveTheory, mode::Symbol; kwargs...)
    if mode == :trace
        contractor = Trace(swt)
        string_formula = "Tr S"
    elseif mode == :perp
        contractor = DipoleFactor(swt)
        string_formula = "∑_ij (I - Q⊗Q){i,j} S{i,j}\n\n(i,j = Sx,Sy,Sz)"
    elseif mode == :full
        contractor = FullTensor(swt)
        string_formula = "S{α,β}"
    end
    intensity_formula(swt,contractor;string_formula,kwargs...)
end

function intensity_formula(swt::SpinWaveTheory, contractor::Contraction; kwargs...)
    return_type = contraction_return_type(contractor)
    intensity_formula(swt,required_correlations(contractor); return_type = return_type,kwargs...) do k,ω,correlations
        intensity = contract(correlations, k, contractor)
    end
end

function intensity_formula(f::Function,swt::SpinWaveTheory,corr_ix::AbstractVector{Int64}; kT = Inf, formfactors = nothing, return_type = Float64, string_formula = "f(Q,ω,S{α,β}[ix_q,ix_ω])")
    (; sys, positions_chem, s̃_mat) = swt
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N  = Nf + 1
    nmodes  = Nf * Nm 
    M = sys.mode == :SUN ? 1 : (Ns-1) # scaling factor (=1) if in the fundamental representation
    sqrt_M = √M
    sqrt_Nm_inv = 1.0 / √Nm

    # Preallocation
    Hmat = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Vmat = zeros(ComplexF64, 2*nmodes, 2*nmodes)
    Avec_pref = zeros(ComplexF64, Nm)
    disp = zeros(Float64, nmodes)
    intensity = zeros(Float64, nmodes)

    # Calculate DSSF 
    formula = function(swt::SpinWaveTheory,q::Vec3)
        _, qmag = chemical_to_magnetic(swt, q)

        swt_hamiltonian!(swt, qmag, Hmat)
        bogoliubov!(disp, Vmat, Hmat, swt.energy_tol)

        for site = 1:Nm
            # note that d is the chemical coordinates
            chemical_coor = positions_chem[site]
            phase = exp(-2im * π  * dot(q, chemical_coor))
            Avec_pref[site] = sqrt_Nm_inv * phase * sqrt_M
        end

        for band = 1:nmodes
            v = Vmat[:, band]
            Avec = zeros(ComplexF64, 3)
            for site = 1:Nm
                @views tS_μ = s̃_mat[:, :, :, site]
                for μ = 1:3
                    for α = 2:N
                        Avec[μ] += Avec_pref[site] * (tS_μ[α, 1, μ] * v[(site-1)*(N-1)+α-1+nmodes] + tS_μ[1, α, μ] * v[(site-1)*(N-1)+α-1])
                    end
                end
            end

            # DD: Generalize this based on list of arbitrary operators, optimize out symmetry, etc.
            Sαβ = Matrix{ComplexF64}(undef,3,3)
            Sαβ[1,1] = real(Avec[1] * conj(Avec[1]))
            Sαβ[1,2] = Avec[1] * conj(Avec[2])
            Sαβ[1,3] = Avec[1] * conj(Avec[3])
            Sαβ[2,2] = real(Avec[2] * conj(Avec[2]))
            Sαβ[2,3] = Avec[2] * conj(Avec[3])
            Sαβ[3,3] = real(Avec[3] * conj(Avec[3]))
            Sαβ[2,1] = conj(Sαβ[1,2]) 
            Sαβ[3,1] = conj(Sαβ[3,1]) 
            Sαβ[3,2] = conj(Sαβ[2,3])

            k = swt.recipvecs_chem * q
            intensity[band] = f(k,disp[band],Sαβ[corr_ix])
        end
        return disp, intensity
    end
    SpinWaveIntensityFormula{return_type}(kT,formfactors,string_formula,formula)
end


