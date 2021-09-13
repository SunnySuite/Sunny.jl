import FastDipole: Vec3, Mat3
import FastDipole.Symmetry: SymOp, Crystal, Bond, canonical_bonds, distance, subcrystal, print_bond_table, lattice_vectors
import FastDipole.Symmetry as S
using LinearAlgebra


### Utility functions

function clean_digits(x, n_digits)
    # keep only n_digits past the decimal point
    x = 10.0^-n_digits * round(10.0^n_digits * x)
    # map -0 to 0
    x == -0.0 ? 0.0 : x
end

# Count the number of symmetries with even/odd parity
function count_symmetries(sym)
    n1 = length(filter(x ->  x[2], sym))
    n2 = length(filter(x -> !x[2], sym))
    return n1, n2
end


### Test construction of diamond lattice

# Spglib inferred symmetry
lat_vecs = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(1, 1, 1), Vec3(-1, -1, -1)] / 8
species = ["C", "C"]
cryst = Crystal(lat_vecs, positions, species)
cbonds = canonical_bonds(cryst, 2.)
dist1 = [distance(cryst, b) for b=cbonds]

# Using explicit symops
base_positions = [Vec3(1, 1, 1) / 8]
base_species = ["C"]
cryst = Crystal(lat_vecs, base_positions, base_species, cryst.symops)
cbonds = canonical_bonds(cryst, 2.)
dist2 = [distance(cryst, b) for b=cbonds]

# Using Hall number
lat_vecs = Mat3(I) # must switch to standard cubic unit cell
base_positions = [Vec3(1, 1, 1) / 4]
@assert cryst.hall_number == 525
cryst = Crystal(lat_vecs, base_positions, base_species, 525)
cbonds = canonical_bonds(cryst, 2.)
dist3 = [distance(cryst, b) for b=cbonds]

# Using international symbol
cryst = Crystal(lat_vecs, base_positions, base_species, "F d -3 m")
cbonds = canonical_bonds(cryst, 2.)
dist4 = [distance(cryst, b) for b=cbonds]

@assert dist1 ≈ dist2 ≈ dist3 ≈ dist4



### FCC lattice, primitive unit cell

lat_vecs = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(0., 0, 0)]
species = ["A"]
cryst = Crystal(lat_vecs, positions, species)
print_bond_table(cryst, 2.)

# Calculate interaction table
cbonds = canonical_bonds(cryst, 2.)
b = cbonds[2]
basis = S.basis_for_symmetry_allowed_couplings(cryst, b)
J = basis' * randn(length(basis))
(bs, Js) = S.all_symmetry_related_interactions(cryst, b, J)
@assert length(Js) == S.bond_multiplicity(cryst, b)


### Triangular lattice, primitive unit cell

lat_vecs = Mat3(1, 0, 0,   1/2, √3/2, 0,   0, 0, 10)
positions = [Vec3(0., 0, 0)]
species = [1]
cryst = Crystal(lat_vecs, positions, species)

print_bond_table(cryst, 5.)


### Kagome lattice

lat_vecs = Mat3(1, 0, 0,   1/2, √3/2, 0,   0, 0, 10)
positions = [Vec3(0, 0, 0), Vec3(0.5, 0, 0), Vec3(0, 0.5, 0)]
species = ["A", "A", "A"]
cryst = Crystal(lat_vecs, positions, species)

print_bond_table(cryst, 3.)


### Arbitrary monoclinic

lat_vecs = lattice_vectors(6, 7, 8, 90, 90, 40)
basis_atoms = [Vec3(0,0,0)]
basis_labels = ["A"]
cryst = Crystal(lat_vecs,basis_atoms,basis_labels,"C 2/c")
display(cryst)


### Arbitrary trigonal

lat_vecs = lattice_vectors(5, 5, 6, 90, 90, 120)
basis_atoms = [Vec3(0,0,0)]
basis_labels = ["A"]
cryst = Crystal(lat_vecs,basis_atoms,basis_labels,"P -3")
cryst = Crystal(lat_vecs,basis_atoms,basis_labels,"I:147") # international number
cryst = Crystal(lat_vecs,basis_atoms,basis_labels,"147") # international number
cryst = Crystal(lat_vecs,basis_atoms,basis_labels, 435) # Hall number
cryst = Crystal(lat_vecs,basis_atoms,basis_labels,"R -3")
display(cryst)


### Arbitrary triclinic

lat_vecs = lattice_vectors(6, 7, 8, 70, 80, 90)
basis_atoms = [Vec3(0,0,0)]
basis_labels = ["A"]
cryst = Crystal(lat_vecs, basis_atoms, basis_labels, "P 1")
cryst = Crystal(lat_vecs, basis_atoms, basis_labels) # Infers 'P -1'
display(cryst)
print_bond_table(cryst, 8.)


### Print FeI2 bond tables

using StaticArrays

cryst = subcrystal(Crystal("/Users/kbarros/Desktop/cifs/FeI2.cif"), "Fe2+")
display(cryst)
print_bond_table(cryst, 6.)

cryst = Crystal("/Users/kbarros/Desktop/cifs/diamond_Nature1958.cif")
print_bond_table(cryst, 5.)
