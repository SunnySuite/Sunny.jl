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

Base.@deprecate suggest_magnetic_supercell(qs, dims) suggest_magnetic_supercell(qs)
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


Base.@deprecate lorentzian06(x, η) let
    @warn "`lorentzian06(x, η)` is deprecated! Use `lorentzian06(; fwhm=2η)(x)` instead."
    return lorentzian06(; fwhm=2η)(x)
end

Base.@deprecate lorentzian06(η) let
    @warn "`lorentzian06(η)` is deprecated! Use `lorentzian06(; fwhm=2η)` instead."
    return lorentzian06(; fwhm=2η)
end

Base.@deprecate integrated_lorentzian(η::Float64) let
    @warn "`integrated_lorentzian(η)` is deprecated! Use `integrated_lorentzian(; fwhm=2η)` instead."
    return integrated_lorentzian(; fwhm=2η)
end

Base.@deprecate set_external_field!(sys::System, B) let
    @warn "`set_external_field!(sys, B)` is deprecated! Consider `set_field!(sys, B*units.T)` where `units = Units(:meV, :angstrom)`."
    set_field!(sys, Vec3(B) * Units(:meV, :angstrom).T)
end

Base.@deprecate set_external_field_at!(sys::System, B, site) let
    @warn "`set_external_field_at!(sys, B, site)` is deprecated! Consider `set_field_at!(sys, B*units.T, site)` where `units = Units(:meV, :angstrom)`."
    set_field_at!(sys, Vec3(B) * Units(:meV, :angstrom).T, site)
end

Base.@deprecate rotation_in_rlu(cryst::Crystal, axis, angle) let
    @warn "`rotation_in_rlu(sys, axis, angle)` is deprecated! Consider `domain_average` instead."
    return rotation_in_rlu(cryst, (axis, angle))
end

# Consider `units.K` where `units = Units(:meV, :angstrom)`.
Base.@deprecate_binding meV_per_K Units(:meV, :angstrom).K


function Base.getproperty(x::Type{Units}, name::Symbol)
    if name in (:theory, :meV)
        @warn "Units.$name is deprecated! See `Units` docs for new interface."
        return nothing
    end
    return getfield(x, name)
end

Base.@deprecate dynamic_correlations(sys; opts...) let
    error("Use SampledCorrelations(...) instead of dynamic_correlations(...)")
end

Base.@deprecate instant_correlations(sys; opts...) let
    error("Use SampledCorrelationsStatic(...) instead of instant_correlations(...)")
end

Base.@deprecate intensity_formula(opts1...; opts2...) let
    error("Formulas have been removed. Call `intensities` directly.")
end

Base.@deprecate reciprocal_space_path(cryst::Crystal, qs, density) let
    error("Use q_space_path(...) instead of reciprocal_space_path(...)")
end

Base.@deprecate set_spiral_order_on_sublattice!(sys, i; q=nothing, k=nothing, axis, S0) let
    error("Use repeat_periodically_as_spiral(...) instead of set_spiral_order_on_sublattice!(...)")
end

Base.@deprecate set_spiral_order!(sys; q=nothing, k=nothing, axis, S0) let
    error("Use repeat_periodically_as_spiral(...) instead of set_spiral_order!(...)")
end

Base.@deprecate System(crystal::Crystal, dims::NTuple{3,Int}, infos, mode::Symbol; seed=nothing, units=nothing) let
    @warn "Deprecation warning! `dims` is now a keyword argument, e.g., System(...; dims=$dims)"
    return System(crystal, infos, mode; dims, seed, units)
end

Base.@deprecate SpinInfo(i; S, g) let
    @warn "SpinInfo(i; S, g) is deprecated! Use `i => Moment(; S, g)` instead."
    i => Moment(; S, g)
end

# REMEMBER TO ALSO DELETE:
#
# * view_crystal(cryst, max_dist)
# * dims argument in view_crystal and plot_spins
# * λ argument in Langevin constructor
# * Δt argument in dynamic_correlations
# * large_S argument in set_exchange! and set_exchange_at!
# * Argument q in set_spiral_order
# * Argument units to System
# * Missing μ0_μB² in enable_dipole_dipole! and
#   modify_exchange_with_truncated_dipole_dipole!
# * energy_ϵ argument in SpinWaveTheory
