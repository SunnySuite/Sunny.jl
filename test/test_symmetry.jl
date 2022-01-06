@testset "Crystal Construction" begin

### Test construction of diamond lattice

# Spglib inferred symmetry
lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
positions = [[1, 1, 1], [-1, -1, -1]] / 8
cryst = Crystal(lat_vecs, positions)
ref_bonds = reference_bonds(cryst, 2.)
dist1 = [distance(cryst, b) for b in ref_bonds]

# Using explicit symops
lat_vecs = Sunny.Mat3(lat_vecs)
positions = [Sunny.Vec3(1, 1, 1) / 8]
types = [""]
cryst = Sunny.crystal_from_symops(lat_vecs, positions, types, cryst.symops, cryst.spacegroup)
ref_bonds = reference_bonds(cryst, 2.)
dist2 = [distance(cryst, b) for b in ref_bonds]

# Using Hall number
lat_vecs = lattice_vectors(1, 1, 1, 90, 90, 90) # must switch to standard cubic unit cell
positions = [Sunny.Vec3(1, 1, 1) / 4]
cryst = Sunny.crystal_from_hall_number(lat_vecs, positions, types, 525)
ref_bonds = reference_bonds(cryst, 2.)
dist3 = [distance(cryst, b) for b in ref_bonds]

# Using international symbol
positions = [[1, 1, 1] / 4]
# cryst = Crystal(lat_vecs, positions, "F d -3 m") # Ambiguous!
cryst = Crystal(lat_vecs, positions, "F d -3 m"; setting="1")
ref_bonds = reference_bonds(cryst, 2.)
dist4 = [distance(cryst, b) for b in ref_bonds]

@test dist1 ≈ dist2 ≈ dist3 ≈ dist4



### FCC lattice, primitive vs. standard unit cell

lat_vecs = [1 1 0; 1 0 1; 0 1 1]' / 2
positions = [[0, 0, 0]]
cryst = Crystal(lat_vecs, positions)

lat_vecs = [1 0 0; 0 1 0; 0 0 1]'
positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
cryst′ = Crystal(lat_vecs, positions)

@test cryst.sitesyms[1] == cryst′.sitesyms[1]

# Calculate interaction table
ref_bonds = reference_bonds(cryst, 2.)
b = ref_bonds[2]
basis = basis_for_symmetry_allowed_couplings(cryst, b)
J = basis' * randn(length(basis))
(bs, Js) = all_symmetry_related_couplings_for_atom(cryst, b.i, b, J)
@test length(Js) == coordination_number(cryst, b.i, b)


### Triangular lattice, primitive unit cell

c = 10
lat_vecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
positions = [[0, 0, 0]]
cryst = Crystal(lat_vecs, positions)
@test cell_type(cryst) == Sunny.hexagonal
@test nbasis(cryst) == 1
@test cell_volume(cryst) ≈ c * √3 / 2 
@test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))

### Kagome lattice

lat_vecs = [1 0 0;  -1/2 √3/2 0;  0 0 c]'
positions = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
cryst = Crystal(lat_vecs, positions)
@test cell_type(cryst) == Sunny.hexagonal
@test nbasis(cryst) == 3
@test cell_volume(cryst) ≈ c * √3 / 2 
@test all(lattice_params(cryst) .≈ (1., 1., c, 90., 90., 120.))


### Arbitrary monoclinic

mono_lat_params = (6, 7, 8, 90, 90, 40)
lat_vecs = lattice_vectors(mono_lat_params...)
positions = [[0,0,0]]
# cryst = Crystal(lat_vecs, positions, "C 2/c")
cryst = Crystal(lat_vecs, positions, "C 2/c", setting="c1")
@test cell_type(cryst) == Sunny.monoclinic
@test nbasis(cryst) == 4
@test all(lattice_params(cryst) .≈ mono_lat_params)


### Arbitrary trigonal

lat_vecs = lattice_vectors(5, 5, 6, 90, 90, 120)
positions = [[0,0,0]]
cryst1 = Crystal(lat_vecs, positions, "P -3")
@test nbasis(cryst1) == 1
@test cell_type(cryst1) == Sunny.hexagonal
cryst2 = Crystal(lat_vecs, positions, "R -3")
@test nbasis(cryst2) == 3
cryst3 = Crystal(lat_vecs, positions, 147) # spacegroup number
@test cell_type(cryst1) == cell_type(cryst2) == cell_type(cryst3) == Sunny.hexagonal


### Arbitrary triclinic

lat_vecs = lattice_vectors(6, 7, 8, 70, 80, 90)
positions = [[0,0,0]]
cryst1 = Crystal(lat_vecs, positions, "P 1")
@test nbasis(cryst1) == 1
cryst2 = Crystal(lat_vecs, positions) # Infers 'P -1'
@test nbasis(cryst1) == nbasis(cryst2) == 1
@test cell_type(cryst1) == cell_type(cryst2) == Sunny.triclinic

end