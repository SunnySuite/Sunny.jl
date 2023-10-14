Base.@deprecate spin_matrices(; N::Int) let
    @warn "`spin_matrices(; N)` will soon be removed! Use `spin_matrices(spin)` instead where `spin = (N-1)/2`."
    spin_matrices((N-1)/2)
end

Base.@deprecate_binding large_S_spin_operators spin_matrices(Inf)
Base.@deprecate_binding large_S_stevens_operators stevens_matrices(Inf)

Base.@deprecate spin_operators(sys::System, i::Int) let
    @warn "`spin_operators` will soon be removed! Use `spin_matrices(spin)` instead where `spin = spin_irrep_label(sys, i)`."
    spin_matrices(spin_irrep_label(sys, i))
end
Base.@deprecate stevens_operators(sys::System, i::Int) let
    @warn "`stevens_operators` will soon be removed! Use `stevens_matrices(spin)` instead where `spin = spin_irrep_label(sys, i)`."
    stevens_matrices(spin_irrep_label(sys, i))
end

Base.@deprecate suggest_magnetic_supercell(qs, latsize) suggest_magnetic_supercell(qs)
Base.@deprecate offline_viewers() ()
