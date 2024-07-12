Base.@deprecate spin_matrices(; N::Int) let
    @warn "`spin_matrices(; N)` is deprecated! Use `spin_matrices((N-1)/2)` instead."
    spin_matrices((N-1)/2)
end

Base.@deprecate_binding large_S_spin_operators spin_matrices(Inf)
Base.@deprecate_binding large_S_stevens_operators stevens_matrices(Inf)

Base.@deprecate spin_operators(sys::System, i::Int) let
    @warn "`spin_operators` is deprecated! Use `spin_matrices(spin_label(sys, i))` instead."
    spin_matrices(spin_label(sys, i))
end
Base.@deprecate stevens_operators(sys::System, i::Int) let
    @warn "`stevens_operators` is deprecated! Use `stevens_matrices(spin_label(sys, i))` instead."
    stevens_matrices(spin_label(sys, i))
end

Base.@deprecate suggest_magnetic_supercell(qs, latsize) suggest_magnetic_supercell(qs)
Base.@deprecate offline_viewers() ()

function Base.getproperty(value::Langevin, name::Symbol)
    if name == :Δt
        @warn "`Δt` field is deprecated! Use `dt` instead."
        name = :dt
    end
    return getfield(value, name)
end
function Base.setproperty!(value::Langevin, name::Symbol, x)
    if name == :Δt
        @warn "`Δt` field is deprecated! Use `dt` instead."
        name = :dt
    end
    return setfield!(value, name, convert(fieldtype(Langevin, name), x))
end


Base.@deprecate lorentzian(x, η) let
    @warn "`lorentzian(x, η)` is deprecated! Use `lorentzian(; fwhm=2η)(x)` instead."
    return lorentzian(; fwhm=2η)(x)
end

Base.@deprecate lorentzian(η) let
    @warn "`lorentzian(η)` is deprecated! Use `lorentzian(; fwhm=2η)` instead."
    return lorentzian(; fwhm=2η)
end

Base.@deprecate integrated_lorentzian(η::Float64) let
    @warn "`integrated_lorentzian(η)` is deprecated! Use `integrated_lorentzian(; fwhm=2η)` instead."
    return integrated_lorentzian(; fwhm=2η)
end

Base.@deprecate set_external_field!(sys::System, B) let
    @warn "`set_external_field!(sys, B)` is deprecated! Consider `set_field!(sys, B*units.T)` where `units = Units(:meV)`."
    set_field!(sys, Vec3(B) * Units(:meV).T)
end

Base.@deprecate set_external_field_at!(sys::System, B, site) let
    @warn "`set_external_field_at!(sys, B, site)` is deprecated! Consider `set_field_at!(sys, B*units.T, site)` where `units = Units(:meV)`."
    set_field_at!(sys, Vec3(B) * Units(:meV).T, site)
end


# Consider `units.K` where `units = Units(:meV)`.
Base.@deprecate_binding meV_per_K Units(:meV).K


function Base.getproperty(x::Type{Units}, name::Symbol)
    if name in (:theory, :meV)
        @warn "Units.$name is deprecated! See `Units` docs for new interface."
        return nothing
    end
    return getfield(x, name)
end


# REMEMBER TO ALSO DELETE:
#
# * view_crystal(cryst, max_dist)
# * λ argument in Langevin constructor
# * Δt argument in dynamical_correlations
# * large_S argument in set_exchange! and set_exchange_at!
# * Argument q in set_spiral_order*
# * Argument units to System
# * Missing μ0_μB² in enable_dipole_dipole! and
#   modify_exchange_with_truncated_dipole_dipole!
# * hermitianpart[!] for VERSION < v"1.10"
