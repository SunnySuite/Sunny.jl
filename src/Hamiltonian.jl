"""
Defines HamiltonianCPU which maintains the actual internal interaction
types and orchestrates energy/field calculations.
"""

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
    diag_coups  :: Vector{DiagonalCouplingCPU{D}}
    gen_coups   :: Vector{GeneralCouplingCPU{D}}
    dipole_int  :: Union{Nothing, DipoleRealCPU, DipoleFourierCPU}
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
    diag_coups  = Vector{DiagonalCouplingCPU{3}}()
    gen_coups   = Vector{GeneralCouplingCPU{3}}()
    dipole_int  = nothing

    D = dimension(crystal)
    validate_dimensionality(ints, D)

    for int in ints
        if isa(int, ExternalField)
            if isnothing(ext_field)
                ext_field = int
            else
                ext_field.B = ext_field.B + int.B
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
        ext_field, heisenbergs, diag_coups, gen_coups, dipole_int
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
