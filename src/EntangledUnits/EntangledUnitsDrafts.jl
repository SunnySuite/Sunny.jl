################################################################################
# Drafts of alternative approach. Revisit later
################################################################################
# function accum_pair_coupling_into_bond_operator_between_units!(bond_operator, pc, sys, contraction_info)
#     (; bond, scalar, bilin, biquad, general) = pc
#     (; i, j, n) = bond
#     unit1, unitsub1 = contraction_info.forward[i]
#     unit2, unitsub2 = contraction_info.forward[j]
# 
#     Ns_local = Ns_in_units(sys, contraction_info)
#     Ns_contracted = map(Ns -> prod(Ns), Ns_local)
#     Ns1 = Ns_local[unit1]
#     Ns2 = Ns_local[unit2]
#     N1 = sys.Ns[1, 1, 1, i]
#     N2 = sys.Ns[1, 1, 1, j]
#     N = Ns_contracted[unit1] * Ns_contracted[unit2]
#     newbond = Bond(unit1, unit2, n)
# 
#     # bond_operator = zeros(ComplexF64, N, N)
#     @assert size(bond_operator, 1) == size(bond_operator, 2) == N
# 
#     # Add scalar part
#     bond_operator .+= scalar*I(N)
# 
#     # Add bilinear part
#     J = bilin isa Float64 ? bilin*I(3) : bilin
#     Si = [kron(local_op_to_product_space(Sa, unitsub1, Ns1), I(Ns_contracted[unit1])) for Sa in spin_matrices((N1-1)/2)]
#     Sj = [kron(I(Ns_contracted[unit2]), local_op_to_product_space(Sa, unitsub2, Ns2)) for Sa in spin_matrices((N2-1)/2)]
#     bond_operator .+= Si' * J * Sj
# 
#     # Add biquadratic part
#     K = biquad isa Float64 ? diagm(biquad * Sunny.scalar_biquad_metric) : biquad
#     Oi = [kron(local_op_to_product_space(Oa, unitsub1, Ns1), I(Ns_contracted[unit1])) for Oa in stevens_matrices_of_dim(2; N=N1)]
#     Oj = [kron(I(Ns_contracted[unit2]), local_op_to_product_space(Ob, unitsub2, Ns2)) for Ob in stevens_matrices_of_dim(2; N=N2)]
#     bond_operator .+= Oi' * K * Oj
# 
#     # Add general part
#     for (A, B) in general.data
#         bond_operator .+= kron( local_op_to_product_space(A, unitsub1, Ns1), I(Ns_contracted[unit1]) ) * 
#                          kron( I(Ns_contracted[unit2]), local_op_to_product_space(B, unitsub2, Ns2) )
#     end
# 
#     return (; newbond, bond_operator)
# end


# function entangle_system_new(sys::System{M}, units) where M
#     # Assess suitability of system for entangled unit formalism
#     if !isnothing(sys.ewald) 
#         error("Original system has long-range dipole-dipole interactions. These are not currently supported for entangled units.")
#     end
#     if !is_homogeneous(sys)
#         error("Entangled units for inhomogenous systems not yet supported.")
#     end
# 
#     # Construct contracted crystal
#     contracted_crystal, contraction_info = contract_crystal(sys.crystal, units)
# 
#     # Determine Ns for local Hilbert spaces (all must be equal). (TODO: Determine if alternative behavior preferable in mixed case.)
#     Ns_unit = Ns_in_units(sys, contraction_info)
#     Ns_contracted = map(Ns -> prod(Ns), Ns_unit)
#     @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."
# 
#     # Construct empty contracted system
#     dims = size(sys.dipoles)[1:3]
#     spin_infos = [i => Moment(s=(N-1)/2, g=1.0) for (i, N) in enumerate(Ns_contracted)]  # TODO: Decisions about g-factor 
#     sys_entangled = System(contracted_crystal, spin_infos, :SUN; dims)
# 
#     # Construct buffers for all relevant interactions
#     # TODO: Generalize to inhomogenous case (way to avoid excessive allocations for large supercell?)
# 
#     # First determine bonds for which we must create a buffer by mapping the
#     # bonds of original interactions to bonds of entangled system.
#     onsite_operators = [zeros(ComplexF64, Ns_contracted[i], Ns_contracted[i]) for i in 1:natoms(sys_entangled.crystal)]
#     new_bonds = Bond[] 
#     for interaction in sys.interactions_union
#         for pc in interaction.pair
#             (; isculled, bond) = pc
#             if !isculled && !bond_is_in_unit(bond, contraction_info)
#                 (; i, j, n) = bond
#                 i_new, _ = contraction_info.forward[i]
#                 j_new, _ = contraction_info.forward[j]
#                 push!(new_bonds, Bond(i_new, j_new, n))
#             end
#         end
#     end
#     new_bonds = unique(new_bonds)
# 
#     bond_operators = Dict{Bond, Array{ComplexF64, 2}}()
#     for bond in unique(new_bonds)
#         (; i, j) = bond
#         N_bond = Ns_contracted[i]*Ns_contracted[j]
#         bond_operators[bond] = zeros(ComplexF64, N_bond, N_bond)
#     end
# 
# 
#     # Accumulate interactions from original system into buffers
#     for (site_origin, interaction) in enumerate(sys.interactions_union)
#         site_contracted = contraction_info.forward[site_origin][1]
# 
#         # Convert onsite couplings and Zeeman
#         onsite_original = interaction.onsite
#         contracted_site, unit_index = contraction_info.forward[site_origin]
#         Ns = Ns_unit[contracted_site]
#         onsite_operators[contracted_site] += local_op_to_product_space(onsite_original, unit_index, Ns)
# 
#         # Convert pair couplings
#         for pc in interaction.pair
#             pc.isculled && continue
# 
#             bond_original = pc.bond
#             if bond_is_in_unit(bond_original, contraction_info)
#                 accum_pair_coupling_into_bond_operator_in_unit!(onsite_operators[site_contracted], pc, sys, site_contracted, contraction_info)
#             elseif bond_original in keys(bond_operators)
#                 (; i, j, n) = bond_original
#                 new_bond = Bond(contraction_info.forward[i][1], contraction_info.forward[j][1], n)
#                 accum_pair_coupling_into_bond_operator_between_units!(bond_operators[new_bond], pc, sys, contraction_info)
#             end
#         end
#     end
# 
#     # Assign interactions to entangled system
#     for i in 1:natoms(contracted_crystal)
#         set_onsite_coupling!(sys_entangled, onsite_operators[i], i)
#     end
#     
#     for bond in keys(bond_operators)
#         bond_operator = bond_operators[bond]
#         set_pair_coupling!(sys_entangled, bond_operator, bond)
#     end
# 
#     return (; sys_entangled, contraction_info, Ns_unit)
# end