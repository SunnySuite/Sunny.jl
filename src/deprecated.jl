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

function Base.copy(dyn::Langevin)
    @warn "Base.copy(dyn::Langevin) will break in Sunny v0.6! Use `Langevin(dyn.dt; dyn.damping, dyn.kT)` instead."
    Langevin(dyn.dt; dyn.damping, dyn.kT)
end

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


# REMEMBER TO ALSO DELETE:

# view_crystal(cryst, max_dist)
# λ argument in Langevin constructor
# Δt argument in dynamical_correlations
# large_S and biquad arguments in set_exchange! and set_exchange_at!
