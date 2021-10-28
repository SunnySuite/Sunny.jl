# Functions associated with HamiltonianCPU, which maintains the actual internal
# interaction types and orchestrates energy/field calculations.


function validate_and_clean_interactions(ints::Vector{Interaction}, crystal::Crystal, latsize::Vector{Int64})
    D = dimension(crystal)

    # Now that we know dimension D, we can convert every OnSiteQuadratic to
    # QuadraticInteraction
    ints = map(ints) do int
        if isa(int, OnSiteQuadratic)
            return QuadraticInteraction(int.J, Bond{D}(int.site, int.site, zeros(D)), int.label)
        else
            return int
        end
    end

    # Validate all interactions
    for int in ints
        if isa(int, QuadraticInteraction)
            b = int.bond

            # Verify that the dimension is correct
            if length(b.n) != D
                error("Interaction $(repr(MIME("text/plain"), int)) inconsistent with crystal dimension $D.")
            end

            # Verify that the interactions are symmetry-consistent
            if !is_coupling_valid(crystal, b, int.J)
                println("Symmetry-violating interaction: $(repr(MIME("text/plain"), int)).")
                if b.i == b.j && iszero(b.n)
                    println("Allowed single-ion anisotropy for this atom:")
                else
                    println("Allowed exchange for this bond:")
                end
                print_allowed_coupling(crystal, b; prefix="    ")
                println("Use `print_bond(crystal, bond)` for more information.")
                error("Interaction violates symmetry.")
            end

            # Verify that no bond wraps the entire system
            bs = all_symmetry_related_bonds(crystal, b)
            wraps = any(bs) do b
                any(abs.(b.n) .>= latsize)
            end
            if wraps
                println("Distance-violating interaction: $int.")
                error("Interaction wraps system.")
            end

        elseif isa(int, DipoleDipole)
            if D != 3
                error("Dipole-dipole interactions require three dimensions.")
            end
        end
    end

    return ints
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
function HamiltonianCPU(ints::Vector{<:Interaction}, crystal::Crystal, latsize::Vector{Int64})
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCPU{3}}()
    diag_coups  = Vector{DiagonalCouplingCPU{3}}()
    gen_coups   = Vector{GeneralCouplingCPU{3}}()
    dipole_int  = nothing

    ints = validate_and_clean_interactions(ints, crystal, latsize)

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
