"""
    HamiltonianCUDA{D}

Stores and orchestrates the types that perform the actual implementations
of all interactions internally. CUDA Version.
"""
struct HamiltonianCUDA{D}
    ext_field   :: Union{Nothing, ExternalFieldCUDA}
    heisenbergs :: Vector{HeisenbergCUDA{D}}
    diag_coups  :: Vector{DiagonalCouplingCUDA{D}}
    gen_coups   :: Vector{GeneralCouplingCUDA{D}}
    dipole_int  :: Union{Nothing, DipoleFourierCUDA}
    spin_mags   :: CUDA.CuVector{Float64}
end

"""
    HamiltonianCUDA(ints::Vector{<:Interaction}, crystal, latsize, sites_info::Vector{SiteInfo})

Construct a `HamiltonianCUDA{3}` from a list of interactions, converting
each of the interactions into the proper backend type specialized
for the given `crystal` and `latsize`.

Note that `sites_info` must be complete when passed to this constructor.
"""
function HamiltonianCUDA(ints::Vector{<:AbstractInteraction}, crystal::Crystal, latsize::Vector{Int64}, sites_info::Vector{SiteInfo})
    ext_field   = nothing
    heisenbergs = Vector{HeisenbergCUDA{3}}()
    diag_coups  = Vector{DiagonalCouplingCUDA{3}}()
    gen_coups   = Vector{GeneralCouplingCUDA{3}}()
    dipole_int  = nothing
    spin_mags = [site.S for site in sites_info]

    ints = validate_and_clean_interactions(ints, crystal, latsize)

    for int in ints
        if isa(int, ExternalField)
            if isnothing(ext_field)
                ext_field = ExternalFieldCUDA(int, sites_info)
            else
                ext_field.Bgs .+= ExternalFieldCUDA(int, sites_info).Bgs
            end
        elseif isa(int, QuadraticInteraction)
            int_impl = convert_quadratic_cuda(int, crystal, sites_info)
            if isa(int_impl, HeisenbergCUDA)
                push!(heisenbergs, int_impl)
            elseif isa(int_impl, DiagonalCouplingCUDA)
                push!(diag_coups, int_impl)
            elseif isa(int_impl, GeneralCouplingCUDA)
                push!(gen_coups, int_impl)
            else
                error("Quadratic interaction failed to convert to known backend type.")
            end
        elseif isa(int, DipoleDipole)
            if !isnothing(dipole_int)
                @warn "Provided multiple dipole interactions. Only using last one."
            end
            dipole_int = DipoleFourierCUDA(int, crystal, latsize, sites_info)
        else
            error("$(int) failed to convert to known backend type.")
        end
    end

    return HamiltonianCUDA{3}(
        ext_field, heisenbergs, diag_coups, gen_coups, dipole_int, spin_mags
    )
end

# Note: These functions only work for HamiltonianCUDA{3}, since kernels are specialized
#        to 3D.

function energy(spins::CUDA.CuArray{Vec3}, ℋ::HamiltonianCUDA{3}) :: Float64
    E = 0.0
    # Each one of these incurs a kernel launch, but they can be done asynchronously.
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

"""
Updates `B` in-place to hold the local field on `spins` under `ℋ`,
defined as:

``𝐁_i = -∇_{𝐬_i} ℋ / S_i``

with ``𝐬_i`` the unit-vector variable at site i, and ``S_i`` is
the magnitude of the associated spin.

Note that all `_accum_neggrad!` functions should return _just_ the
``-∇_{𝐬_i} ℋ`` term, as the scaling by spin magnitude happens in
this function. Likewise, all code which utilizes local fields should
be calling _this_ function, not the `_accum_neggrad!`'s directly.
"""
function field!(B::CUDA.CuArray{Vec3}, spins::CUDA.CuArray{Vec3}, ℋ::HamiltonianCUDA{3})
    fill!(B, SA[0.0, 0.0, 0.0])
    # Each one of these incurs a kernel launch, but they can be done asynchronously.
    if !isnothing(ℋ.ext_field)
        _accum_neggrad!(B, ℋ.ext_field)
    end
    for heisen in ℋ.heisenbergs
        _accum_neggrad!(B, spins, heisen)
    end
    for diag_coup in ℋ.diag_coups
        _accum_neggrad!(B, spins, diag_coup)
    end
    for gen_coup in ℋ.gen_coups
        _accum_neggrad!(B, spins, gen_coup)
    end
    if !isnothing(ℋ.dipole_int)
        _accum_neggrad!(B, spins, ℋ.dipole_int)
    end

    # Normalize each gradient by the spin magnitude on that sublattice
    B ./= reshape(ℋ.spin_mags, length(ℋ.spin_mages), ones(Int, 3)...)
end