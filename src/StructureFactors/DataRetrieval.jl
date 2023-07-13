################################################################################
# Basic functions for retrieving 𝒮(𝐪,ω) values
################################################################################

abstract type IntensityFormula end

struct ClassicalIntensityFormula{T} <: IntensityFormula
    kT :: Float64
    formfactors
    string_formula :: String
    calc_intensity :: Function
end

function Base.show(io::IO, formula::ClassicalIntensityFormula{T}) where T
    print(io,"ClassicalIntensityFormula{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", formula::ClassicalIntensityFormula{T}) where T
    printstyled(io, "Classical Scattering Intensity Formula\n";bold=true, color=:underline)

    formula_lines = split(formula.string_formula,'\n')

    intensity_equals = "  Intensity[ix_q,ix_ω] = "
    println(io,"At discrete scattering modes S = S[ix_q,ix_ω], use:")
    println(io)
    println(io,intensity_equals,formula_lines[1])
    for i = 2:length(formula_lines)
        precursor = repeat(' ', textwidth(intensity_equals))
        println(io,precursor,formula_lines[i])
    end
    println(io)

    if isnothing(formula.formfactors)
        printstyled(io, "No form factors specified\n";color=:yellow)
    else
        printstyled(io, "Form factors included in S ✓\n";color=:green)
    end
    if formula.kT == Inf
        printstyled(io, "No temperature correction";color=:yellow)
        print(io, " (kT = ∞)\n")
    else
        printstyled(io, "Temperature corrected (kT = $(formula.kT)) ✓\n";color = :green)
    end
    if T != Float64
        println(io,"Intensity :: $(T)")
    end
end

"""
    formula = intensity_formula(sf::StructureFactor; kwargs...)
    formula.calc_intensity(sf,q,ix_q,ix_ω)

Establish a formula for computing the intensity of the discrete scattering modes `(q,ω)` using the correlation data ``𝒮^{αβ}(q,ω)`` stored in the [`StructureFactor`](@ref).
The `formula` returned from `intensity_formula` can be passed to [`intensities_interpolated`](@ref) or [`intensities_binned`](@ref).

Sunny has several built-in formulas that can be selected by setting `contraction_mode` to one of these values:

- `:perp` (default), which contracts ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} - q_{α}q_{β}``, returning the unpolarized intensity.
- `:trace`, which yields ``\\operatorname{tr} 𝒮(q,ω) = ∑_α 𝒮^{αα}(q,ω)``
- `:full`, which will return all elements ``𝒮^{αβ}(𝐪,ω)`` without contraction.

Additionally, there are keyword arguments providing temperature and form factor corrections:

- `kT`: If a temperature is provided, the intensities will be rescaled by a
    temperature- and ω-dependent classical-to-quantum factor. `kT` should be
    specified when making comparisons with spin wave calculations or
    experimental data. If `kT` is not specified, infinite temperature (no correction) is assumed.
- `formfactors`: To apply form factor corrections, provide this keyword with a
    vector of `FormFactor`s, one for each unique site in the unit cell. The form factors
    will be symmetry propagated to all equivalent sites.

Alternatively, a custom formula can be specifed by providing a function `intensity = f(q,ω,correlations)` and specifying which correlations it requires:

    intensity_formula(f,sf::StructureFactor, required_correlations; kwargs...)

The function is intended to be specified using `do` notation. For example, this custom formula sums the off-diagonal correlations:

    required = [(:Sx,:Sy),(:Sy,:Sz),(:Sx,:Sz)]
    intensity_formula(sf,required,return_type = ComplexF64) do k, ω, off_diagonal_correlations
        sum(off_diagonal_correlations)
    end

If your custom formula returns a type other than `Float64`, use the `return_type` keyword argument to flag this.
"""
function intensity_formula(f::Function,sf::StructureFactor,required_correlations; kwargs...)
    # SQTODO: This corr_ix may contain repeated correlations if the user does a silly
    # thing like [(:Sx,:Sy),(:Sy,:Sx)], and this can technically be optimized so it's
    # not computed twice
    corr_ix = lookup_correlations(sf,required_correlations)
    intensity_formula(f,sf,corr_ix;kwargs...)
end

function intensity_formula(f::Function,sf::StructureFactor,corr_ix::AbstractVector{Int64}; kT = Inf, formfactors = nothing, return_type = Float64, string_formula = "f(Q,ω,S{α,β}[ix_q,ix_ω])")
    # If temperature given, ensure it's greater than 0.0
    if iszero(kT)
        error("`kT` must be greater than zero.")
    end

    ffdata = prepare_form_factors(sf, formfactors)
    NAtoms = size(sf.data)[2]
    NCorr = length(corr_ix)

    ωs_sf = ωs(sf;negative_energies=true)
    formula = function (sf::StructureFactor,k::Vec3,ix_q::CartesianIndex{3},ix_ω::Int64)
        correlations = phase_averaged_elements(view(sf.data,corr_ix,:,:,ix_q,ix_ω), k, sf, ffdata, Val(NCorr), Val(NAtoms))

        ω = ωs_sf[ix_ω]
        intensity = f(k,ω,correlations) * classical_to_quantum(ω, kT)

        # Having this line saves the return_type in the function closure
        # so that it can be read by intensities later
        intensity :: return_type
    end
    ClassicalIntensityFormula{return_type}(kT,formfactors,string_formula,formula)
end

function classical_to_quantum(ω, kT::Float64)
  if kT == Inf
    return 1.0
  end
  if ω > 0
    ω/(kT*(1 - exp(-ω/kT)))
  elseif iszero(ω)
    1.0
  else
    -ω*exp(ω/kT)/(kT*(1 - exp(ω/kT)))
  end
end

function prepare_form_factors(sf, formfactors)
    if isnothing(formfactors)
        cryst = isnothing(sf.origin_crystal) ? sf.crystal : sf.origin_crystal 
        class_indices = [findfirst(==(class_label), cryst.classes) for class_label in unique(cryst.classes)]
        formfactors = [FormFactor{Sunny.EMPTY_FF}(; atom) for atom in class_indices]
    end
    formfactors = upconvert_form_factors(formfactors) # Ensure formfactors have consistent type
    return propagate_form_factors(sf, formfactors)
end


"""
    lorentzian(x, η) 

Returns ``η/(π(x^2 + η^2))``.
"""
lorentzian(x, η) = η/(π*(x^2 + η^2))

"""
    integrated_lorentzian(η) 

Returns ``x \\mapsto atan(x/η)/π`` for use with [`intensities_binned`](@ref).
"""
integrated_lorentzian(η) = x -> atan(x/η)/π

"""
    broaden_energy(sf::StructureFactor, vals, kernel::Function; negative_energies=false)

Performs a real-space convolution along the energy axis of an array of
intensities. Assumes the format of the intensities array corresponds to what
would be returned by [`intensities_interpolated`](@ref). `kernel` must be a function that
takes two numbers: `kernel(ω, ω₀)`, where `ω` is a frequency, and `ω₀` is the
center frequency of the kernel. Sunny provides [`lorentzian`](@ref)
for the most common use case:

```
newvals = broaden_energy(sf, vals, (ω, ω₀) -> lorentzian(ω-ω₀, 0.2))
```
"""
function broaden_energy(sf::StructureFactor, is, kernel::Function; negative_energies=false)
    dims = size(is)
    ωvals = ωs(sf; negative_energies)
    out = zero(is)
    for (ω₀i, ω₀) in enumerate(ωvals)
        for (ωi, ω) in enumerate(ωvals)
            for qi in CartesianIndices(dims[1:end-1])
                out[qi,ωi] += is[qi,ω₀i]*kernel(ω, ω₀)
            end
        end
    end
    return out
end
