"""Defines Hamiltonian, a user-facing Hamiltonian-defining type, as well as HamiltonianCPU
    which maintains the actual internal interaction types and orchestrates energy/field
    calculations.
"""

"""
    Hamiltonian{D}

Defines a Hamiltonian for a `D`-dimensional spin system.
"""
struct Hamiltonian{D}
    interactions :: Vector{<:Interaction}

    # Inner constructor verifies all interactions have compatible dimensionalities
    function Hamiltonian{D}(ints) where {D}
        for int in ints
            if isa(int, PairInt) && typeof(int).parameters[1] != D
                error("One of the provided pair interactions is not the correct dimensionality.")
            elseif isa(int, DipoleDipole) && D != 3
                error("Dipole-dipole interactions only supported for D = 3.")
            end
        end
        new{D}(collect(ints))
    end
end

"""
    Hamiltonian(ints)
    Hamiltonian(ints...)
Constructor for a `Hamiltonian` which attempts to infer dimensionality
from the provided interactions. Will fail if no pair or dipole
interactions are defined.
"""
function Hamiltonian(ints)
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
               Use explicit constructor Hamiltonian{D}.
            """
        )
    else
        Hamiltonian{D}(ints)
    end
end

Hamiltonian(ints::Vararg{<:Interaction}) = Hamiltonian(collect(ints))

"""
Like `Hamiltonian{D}`, but stores and orchestrates the types that perform
the actual implementations of all interactions internally.
"""
struct HamiltonianCPU{D}
    ext_field   :: Union{Nothing, ExternalField}
    heisenbergs :: Vector{HeisenbergCPU{D}}
    on_sites    :: Vector{OnSite}
    diag_coups  :: Vector{DiagonalCouplingCPU{D}}
    gen_coups   :: Vector{GeneralCouplingCPU{D}}
    dipole_int  :: Union{Nothing, DipoleFourierCPU}
end

"""
    HamiltonianCPU{D}(ℋ::Hamiltonian, crystal, lattice)

Construct a `HamiltonianCPU{D}` from a `Hamiltonian{D}`, converting
each of the interactions into the proper backend type specialized
for the given `crystal` and `latsize`.
"""
function HamiltonianCPU{D}(ℋ::Hamiltonian{D}, crystal::Crystal, latsize) where {D}
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCPU{D}}()
    on_sites    = Vector{OnSite}()
    diag_coups  = Vector{DiagonalCouplingCPU{D}}()
    gen_coups   = Vector{GeneralCouplingCPU{D}}()
    dipole_int  = nothing
    for int in ℋ.interactions
        if isa(int, ExternalField)
            if !isnothing(ext_field)
                @warn "Provided multiple external fields. Only using last one."
            end
            ext_field = int
        elseif isa(int, Heisenberg)
            push!(heisenbergs, HeisenbergCPU(int, crystal))
        elseif isa(int, OnSite)
            push!(on_sites, int)
        elseif isa(int, DiagonalCoupling)
            push!(diag_coups, DiagonalCouplingCPU(int, crystal))
        elseif isa(int, GeneralCoupling)
            push!(gen_coups, GeneralCouplingCPU(int, crystal))
        elseif isa(int, DipoleDipole)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = DipoleFourierCPU(int, crystal, latsize)
        end
    end
    return HamiltonianCPU{D}(
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
    for on_site in ℋ.on_sites
        E += energy(spins, on_site)
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
    for on_site in ℋ.on_sites
        _accum_field!(B, spins, on_site)
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
