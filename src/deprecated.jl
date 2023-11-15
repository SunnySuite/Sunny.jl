Base.@deprecate spin_matrices(; N::Int) let
    @warn "`spin_matrices(; N)` will soon be removed! Use `spin_matrices((N-1)/2)` instead."
    spin_matrices((N-1)/2)
end

Base.@deprecate_binding large_S_spin_operators spin_matrices(Inf)
Base.@deprecate_binding large_S_stevens_operators stevens_matrices(Inf)

Base.@deprecate spin_operators(sys::System, i::Int) let
    @warn "`spin_operators` will soon be removed! Use `spin_matrices(spin_label(sys, i))` instead."
    spin_matrices(spin_label(sys, i))
end
Base.@deprecate stevens_operators(sys::System, i::Int) let
    @warn "`stevens_operators` will soon be removed! Use `stevens_matrices(spin_label(sys, i))` instead."
    stevens_matrices(spin_label(sys, i))
end

Base.@deprecate suggest_magnetic_supercell(qs, latsize) suggest_magnetic_supercell(qs)
Base.@deprecate offline_viewers() ()

function Base.copy(dyn::Langevin)
    @warn "Base.copy(dyn::Langevin) will soon be removed! Use `Langevin(dyn.Δt; dyn.λ, dyn.kT)` instead."
    Langevin(dyn.Δt; dyn.λ, dyn.kT)
end
