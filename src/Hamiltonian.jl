"""
Defines HamiltonianCPU which maintains the actual internal interaction
types and orchestrates energy/field calculations.
"""

function infer_dimensionality(ints::Vector{<:Interaction})
    D = nothing
    # Try to infer dimenisonality from interactions
    for int in ints
        if isa(int, PairInt)
            # Sort of hacky -- is there a better way?
            intD = typeof(int).parameters[1]
        elseif isa(int, DipoleDipole)
            intD = 3
        else
            intD = nothing
        end

        if !isnothing(intD)
            if isnothing(D)
                D = intD
            elseif D != intD
                error(
                    """Provided interactions of multiple inconsistent
                        dimensionalities!
                    """
                )
            end
        end
    end

    if isnothing(D)
        error(
            """Could not infer dimensionality from arguments.
               Use explicit constructor HamiltonianCPU{D}.
            """
        )
    end
    D
end

function validate_dimensionality(ints::Vector{<:Interaction}, D::Int)
    for int in ints
        if isa(int, QuadraticInteraction) && typeof(int).parameters[1] != D
            error("One of the provided pair interactions is not the correct dimensionality.")
        elseif isa(int, DipoleDipole) && D != 3
            error("Dipole-dipole interactions only supported for D = 3.")
        end
    end
end

"""
    HamiltonianCPU{D}

Stores and orchestrates the types that perform the actual implementations
of all interactions internally.
"""
struct HamiltonianCPU{D}
    ext_field   :: Union{Nothing, ExternalField}
    heisenbergs :: Vector{HeisenbergCPU{D}}
    on_sites    :: Union{Nothing, OnSiteQuadraticCPU}
    diag_coups  :: Vector{DiagonalCouplingCPU{D}}
    gen_coups   :: Vector{GeneralCouplingCPU{D}}
    dipole_int  :: Union{Nothing, DipoleFourierCPU}
end

"""
    HamiltonianCPU(ints::Vector{<:Interaction}, crystal, latsize)

Construct a `HamiltonianCPU{3}` from a list of interactions, converting
each of the interactions into the proper backend type specialized
for the given `crystal` and `latsize`.
"""
function HamiltonianCPU(ints::Vector{<:Interaction}, crystal::Crystal, latsize)
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCPU{3}}()
    on_sites    = nothing
    diag_coups  = Vector{DiagonalCouplingCPU{3}}()
    gen_coups   = Vector{GeneralCouplingCPU{3}}()
    dipole_int  = nothing

    validate_dimensionality(ints, 3)

    for int in ints
        if isa(int, ExternalField)
            if isnothing(ext_field)
                ext_field = int
            else
                ext_field.B = ext_field.B + int.B
            end
        elseif isa(int, OnSiteQuadratic)
            if isnothing(on_sites)
                on_sites = OnSiteQuadraticCPU(int, crystal)
            else
                @warn "Multiple on-sites detected. Merging, keeping only last label."
                class_on_site = OnSiteQuadraticCPU(int, crystal)
                newJs = [on_sites.Js[i] + class_on_site.Js[i] for i in 1:nbasis(crystal)]
                on_sites = OnSiteQuadraticCPU(newJs, class_on_site.label)
            end
        elseif isa(int, QuadraticInteraction)
            int_impl = convert_quadratic(int, crystal)
            if isa(int_impl, HeisenbergCPU)
                push!(heisenbergs, int_impl)
            elseif isa(int_impl, DiagonalCouplingCPU)
                push!(diag_coups, int_impl)
            elseif isa(int_impl, GeneralCouplingCPU)
                push!(gen_coups, int_impl)
            else
                error("Quadratic interaction failed to convert to known backend type.")
            end
        elseif isa(int, DipoleDipole)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = DipoleFourierCPU(int, crystal, latsize)
        else
            error("$(int) failed to convert to known backend type.")
        end
    end
    return HamiltonianCPU{3}(
        ext_field, heisenbergs, on_sites,
        diag_coups, gen_coups, dipole_int
    )
end

function energy(spins::Array{Vec3}, ℋ::HamiltonianCPU) :: Float64
    E = 0.0
    if !isnothing(ℋ.ext_field)
        E += energy(spins, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        E += energy(spins, heisen)
    end
    if !isnothing(ℋ.on_sites)
        E += energy(spins, ℋ.on_sites)
    end
    for diag_coup in ℋ.diag_coups
        E += energy(spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        E += energy(spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        E += energy(spins, ℋ.dipole_int)
    end
    return E
end

function field!(B::Array{Vec3}, spins::Array{Vec3}, ℋ::HamiltonianCPU)
    fill!(B, SA[0.0, 0.0, 0.0])
    if !isnothing(ℋ.ext_field)
        _accum_field!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        _accum_field!(B, spins, heisen)
    end
    if !isnothing(ℋ.on_sites)
        _accum_field!(B, spins, ℋ.on_sites)
    end
    for diag_coup in ℋ.diag_coups
        _accum_field!(B, spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        _accum_field!(B, spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_field!(B, spins, ℋ.dipole_int)
    end
end
