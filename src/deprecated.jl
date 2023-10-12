Base.@deprecate spin_matrices(; N::Int) let
    @warn "`spin_matrices(; N)` is deprecated, use `spin_matrices(S)` instead."
    spin_matrices((N-1)/2)
end

Base.@deprecate_binding large_S_spin_operators spin_matrices(Inf)
Base.@deprecate_binding large_S_stevens_operators stevens_matrices(Inf)

Base.@deprecate spin_operators(sys::System, i::Int) let
    @warn "`spin_operators` is deprecated, use `spin_matrices` instead."
    spin_matrices(spin_irrep_label(sys, i))
end
Base.@deprecate stevens_operators(sys::System, i::Int) let
    @warn "`stevens_operators` is deprecated, use `stevens_matrices` instead."
    stevens_matrices(spin_irrep_label(sys, i))
end

Base.@deprecate suggest_magnetic_supercell(qs, latsize) suggest_magnetic_supercell(qs)
Base.@deprecate offline_viewers() ()
