################################################################################
# Types
################################################################################
abstract type Contraction{T} end  # T determines type value returned by the contraction 

struct Trace{N} <: Contraction{Float64}
    indices :: SVector{N, Int64}
    unilateral_to_bilateral :: Bool
end

struct DipoleFactor <: Contraction{Float64}
    indices :: SVector{9,Int64}
    # Should we apply the unilateral-to-bilateral S → S + S' symmetrization?
    unilateral_to_bilateral :: Bool
end

struct Element <: Contraction{ComplexF64}
    index :: Int64
end

struct FullTensor{NCorr,NSquare,NObs,NObs2} <: Contraction{SMatrix{NObs, NObs, ComplexF64, NObs2}}
    indices :: SVector{NSquare, Int64}
end

struct AllAvailable{NCorr} <: Contraction{SVector{NCorr, ComplexF64}}
end


################################################################################
# Constructors
################################################################################
function Trace(obs::ObservableInfo;unilateral_to_bilateral = true)
    # Collect all indices for matrix elements 𝒮^αβ where α=β
    indices = Int64[]
    for (ki,i) = obs.observable_ixs
        autocorrelation_index = CartesianIndex(i,i)
        if haskey(obs.correlations,autocorrelation_index)
            push!(indices,obs.correlations[autocorrelation_index])
        else
            problematic_correlation = ki
            error("Can't calculate trace because auto-correlation of the $problematic_correlation observable was not computed.")
        end
    end

    # SQ N.B.: This error doesn't make much sense, does it?
    # So what if they used a different number from the default number of observables?
    # Case in point: If you are doing dipole correlations in SU(N) mode, you're not taking
    # the full trace, and this will error out.

    #=
    total_autocorrelations = N == 0 ? 3 : N*N-1
    if length(indices) != total_autocorrelations
        error("Unexpected number of observables were encounted. Expected $total_autocorrelations but actually have $(length(sc.observables)): $(keys(sc.observable_ixs))")
    end
    =#

    indices = sort(indices)
    Trace(SVector{length(indices), Int64}(indices),unilateral_to_bilateral)
end

function DipoleFactor(obs::ObservableInfo; unilateral_to_bilateral = true, spin_components = [:Sx,:Sy,:Sz])
    # Ensure that the observables themselves are present
    for si in spin_components
        if !haskey(obs.observable_ixs,si)
            error("Observable $(si) missing, but required for dipole correction factor")
        end
    end

    # Ensure that the required correlations are also present
    sx,sy,sz = spin_components
    dipole_correlations = [(sx,sx),(sy,sx),(sz,sx),(sx,sy),(sy,sy),(sz,sy),(sx,sz),(sy,sz),(sz,sz)]
    indices = lookup_correlations(obs,dipole_correlations; err_msg = αβ -> "Missing correlation $(αβ), which is required to compute the depolarization correction.")
    DipoleFactor(indices, unilateral_to_bilateral)
end

function Element(obs::ObservableInfo, pair::Tuple{Symbol,Symbol})
    Element(only(lookup_correlations(obs,[pair]; err_msg = pair -> "Missing correlation $(pair), which was requested.")))
end

function FullTensor(obs::ObservableInfo)
    n_obs = num_observables(obs)
    tensor_elements = Matrix{Tuple{Symbol,Symbol}}(undef,n_obs,n_obs)
    for (ki,i) = obs.observable_ixs, (kj,j) = obs.observable_ixs
      tensor_elements[i,j] = (ki,kj) # Required to put matrix in correct order
    end
    indices = lookup_correlations(obs, collect(tensor_elements); err_msg = αβ -> "Missing correlation $(αβ). All correlations are required to return the full tensor.")
    FullTensor{num_correlations(obs),length(indices),n_obs,n_obs*n_obs}(indices)
end

################################################################################
# Contraction helper functions
################################################################################
@inline function polarization_matrix(k::Vec3)
    k /= norm(k) + 1e-12
    return SMatrix{3, 3, Float64, 9}(I(3) - k * k')
end

################################################################################
# Contraction methods
#
# The first argument to `contract' is the correlations which were requested
# by required_correlations. In particular, for DipoleFactor and FullTensor,
# these correlations are already sorted in the desired order, so all that needs to
# be done is the actual calculation (no re-ordering)
################################################################################


# Diagonal elements should be real only. Finite imaginary component is 
# usually on order 1e-17 and is due to roundoff in phase_averaged_elements.
function contract(diagonal_elements, _, traceinfo::Trace)
    if traceinfo.unilateral_to_bilateral
        2 * sum(real(diagonal_elements))
    else
        sum(real(diagonal_elements))
    end
end

function contract(dipole_elements, k::Vec3, dipoleinfo::DipoleFactor)
    dip_factor = polarization_matrix(k)

    # Note, can just take the real part since:
    #   (1) diagonal elements are real by construction, and 
    #   (2) pairs of off diagonal contributions have the form x*conj(y) + conj(x)*y = 2real(x*conj(y)).

    Sab = reshape(dipole_elements,3,3)

    #display(Sab)
    #println("[CS] part:")
    #display(Sab + Sab')
    #println("[-CS] part:")
    #display(Sab - Sab')

    # If Sab is the *unilateral (Laplace) transform* of the time-domin
    # correlations, then we are *not* gaurunteed conjugate-symmetry of the matrix,
    # so we need to explicitly symmterize it now to ensure that the result
    # of contracting with the dipole factor is real. This symmetrization corresponds
    # to glueing together the two halves of the unilateral transform to recover
    # the *bilateral (Fourier) transform* of the time-domain correlations.
    # Because of this glueing, there is no factor 1/2 needed.
    #
    # Warning: This procedure assumes real observables, A = A†, since
    # the correct symmterization is over simultaneous
    #   (1) complex conjugation
    #   (2) transpose (= swap observables)
    #   (3) dagger each observable
    # but we only perform steps (1) and (2).
    if dipoleinfo.unilateral_to_bilateral
        Sab = Sab + Sab'
    end

    # Since dip_factor is symmteric, and Sab is (now) guaranteed
    # to be conjugate-symmetric, dipole_intensity is real up to
    # machine precision
    dipole_intensity = sum(dip_factor .* Sab)

    # This assertation catches the case where the user set
    # `unilateral_to_bilateral = false' even though their data
    # is unilateral. In this case, their Sab is not conjugate-symmetric,
    # and we didn't conjugate-symmetrize it, so dipole_intensity is
    # errantly complex
    @assert abs(imag(dipole_intensity)) < 1e-12

    return real(dipole_intensity)
end


contract(specific_element, _, ::Element) = only(specific_element)

function contract(all_elems, _, full::FullTensor{NCorr,NSquare,NObs,NObs2}) where {NCorr, NSquare,NObs,NObs2}
    reshape(all_elems,NObs,NObs)
end

function contract(all_elems, _, ::AllAvailable{NCorr}) where NCorr
    all_elems
end

# The required_correlations specifies which correlations
# (as numbered according to the *values* of the ObservableInfo.correlations map)
# are needed for the calculation, and what order they should be retrieved in.
# This is particularly important for the DipoleFactor and FullTensor,
# which retreive the correlations in a specific order.
required_correlations(traceinfo::Trace) = traceinfo.indices
required_correlations(dipoleinfo::DipoleFactor) = dipoleinfo.indices
required_correlations(eleminfo::Element) = [eleminfo.index]
required_correlations(fullinfo::FullTensor{NCorr}) where NCorr = fullinfo.indices
required_correlations(::AllAvailable{NCorr}) where NCorr = 1:NCorr


################################################################################
# Contraction utils
################################################################################
Base.zeros(::Contraction{T}, dims...) where T = zeros(T, dims...)

function contractor_from_mode(source, mode::Symbol)
    if mode == :trace
        contractor = Trace(source.observables; unilateral_to_bilateral = true)
        string_formula = "Tr S′\n\n with S′ = S + S′"
    elseif mode == :perp
        contractor = DipoleFactor(source.observables; unilateral_to_bilateral = true)
        string_formula = "∑_ij (I - Q⊗Q){i,j} S′{i,j}\n\n(i,j = Sx,Sy,Sz) and with S′ = S + S†"
    elseif mode == :full
        contractor = FullTensor(source.observables)
        string_formula = "S{α,β}"
    elseif mode == :all_available
        corrs = keys(source.observables.correlations)
        contractor = AllAvailable{length(corrs)}()
        string_formula = "[" * join(map(x -> "S{$(x.I[1]),$(x.I[2])}",corrs),", ") * "]"
    else
        error("Unknown mode: $mode. Try one of :trace, :perp, :full, :all_available")
    end
    return contractor, string_formula
end

"""
    intensity_formula([swt or sc], contraction_mode::Symbol)

Sunny has several built-in formulas that can be selected by setting `contraction_mode` to one of these values:

- `:trace` (default), which yields ``\\operatorname{tr} 𝒮(q,ω) = ∑_α 𝒮^{αα}(q,ω)``
- `:perp`, which contracts ``𝒮^{αβ}(q,ω)`` with the dipole factor ``δ_{αβ} - q_{α}q_{β}``, returning the unpolarized intensity.
- `:full`, which will return all elements ``𝒮^{αβ}(𝐪,ω)`` without contraction.
- `:all_available`, which will return each computed correlation as determined by `sc.observables` or `swt.observables`.
"""
function intensity_formula(swt::SpinWaveTheory, mode::Symbol; kwargs...)
    contractor, string_formula = contractor_from_mode(swt, mode)
    intensity_formula(swt, contractor; string_formula, kwargs...)
end

function intensity_formula(swt::SpinWaveTheory, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(swt,required_correlations(contractor); return_type = T,kwargs...) do k,ω,correlations
        intensity = contract(correlations, k, contractor)
    end
end

function intensity_formula(sc::SampledCorrelations, elem::Tuple{Symbol,Symbol}; kwargs...)
    string_formula = "S{$(elem[1]),$(elem[2])}[ix_q,ix_ω]"
    intensity_formula(sc,Element(sc, elem); string_formula, kwargs...)
end

function intensity_formula(sc::SampledCorrelations, mode::Symbol; kwargs...)
    contractor, string_formula = contractor_from_mode(sc, mode)
    intensity_formula(sc, contractor; string_formula, kwargs...)
end

function intensity_formula(sc::SampledCorrelations, contractor::Contraction{T}; kwargs...) where T
    intensity_formula(sc,required_correlations(contractor); return_type = T,kwargs...) do k,ω,correlations
        intensity = contract(correlations, k, contractor)
    end
end
