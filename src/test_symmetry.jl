import FastDipole: Vec3, Mat3
import FastDipole.Symmetry: SymOp, Cell, Bond, canonical_bonds, distance
import FastDipole.Symmetry as S

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


### Diamond lattice

# Spglib inferred symmetry
lattice = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(1, 1, 1), Vec3(-1, -1, -1)] / 8
species = ["C", "C"]
cell = Cell(lattice, positions, species)
cbonds = canonical_bonds(cell, 2.)
dist1 = [distance(cell, b) for b=cbonds]

# explicit symops
base_positions = [Vec3(1, 1, 1) / 8]
base_species = ["C"]
symops = SymOp[SymOp([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [0.0, 0.0, -0.0]), SymOp([0.0 0.0 -1.0; -1.0 0.0 0.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([-1.0 -1.0 -1.0; 0.0 0.0 1.0; 0.0 1.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 -1.0 0.0; 1.0 1.0 1.0; -1.0 0.0 0.0], [0.0, 0.5, -0.0]), SymOp([0.0 1.0 0.0; 1.0 0.0 0.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([-1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 -1.0 0.0], [0.0, 0.0, -0.0]), SymOp([0.0 0.0 1.0; -1.0 -1.0 -1.0; 1.0 0.0 0.0], [0.0, 0.5, -0.0]), SymOp([1.0 1.0 1.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], [0.5, 0.0, -0.0]), SymOp([0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0], [0.0, 0.0, -0.0]), SymOp([-1.0 0.0 0.0; 1.0 1.0 1.0; 0.0 0.0 -1.0], [0.0, 0.5, -0.0]), SymOp([0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([1.0 1.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0], [0.5, 0.0, -0.0]), SymOp([1.0 0.0 0.0; -1.0 -1.0 -1.0; 0.0 1.0 0.0], [0.0, 0.5, -0.0]), SymOp([0.0 0.0 -1.0; 0.0 -1.0 0.0; -1.0 0.0 0.0], [0.0, 0.0, -0.0]), SymOp([-1.0 -1.0 -1.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.5, 0.0, -0.0]), SymOp([0.0 -1.0 0.0; 0.0 0.0 -1.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0], [0.0, 0.0, -0.0]), SymOp([1.0 1.0 1.0; 0.0 0.0 -1.0; -1.0 0.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 1.0 0.0; -1.0 -1.0 -1.0; 0.0 0.0 1.0], [0.0, 0.5, -0.0]), SymOp([-1.0 0.0 0.0; 0.0 -1.0 0.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([-1.0 -1.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 -1.0], [0.0, 0.0, -0.0]), SymOp([1.0 0.0 0.0; 0.0 0.0 1.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([0.0 0.0 -1.0; 1.0 1.0 1.0; 0.0 -1.0 0.0], [0.0, 0.5, -0.0]), SymOp([-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], [0.0, 0.0, -0.0]), SymOp([0.0 0.0 1.0; 1.0 0.0 0.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([1.0 1.0 1.0; 0.0 0.0 -1.0; 0.0 -1.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 1.0 0.0; -1.0 -1.0 -1.0; 1.0 0.0 0.0], [0.0, 0.5, -0.0]), SymOp([0.0 -1.0 0.0; -1.0 0.0 0.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0], [0.0, 0.0, -0.0]), SymOp([0.0 0.0 -1.0; 1.0 1.0 1.0; -1.0 0.0 0.0], [0.0, 0.5, -0.0]), SymOp([-1.0 -1.0 -1.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [0.5, 0.0, -0.0]), SymOp([0.0 -1.0 0.0; 0.0 0.0 -1.0; -1.0 0.0 0.0], [0.0, 0.0, -0.0]), SymOp([1.0 0.0 0.0; -1.0 -1.0 -1.0; 0.0 0.0 1.0], [0.0, 0.5, -0.0]), SymOp([0.0 0.0 -1.0; 0.0 -1.0 0.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([-1.0 -1.0 -1.0; 1.0 0.0 0.0; 0.0 1.0 0.0], [0.5, 0.0, -0.0]), SymOp([-1.0 0.0 0.0; 1.0 1.0 1.0; 0.0 -1.0 0.0], [0.0, 0.5, -0.0]), SymOp([0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0], [0.0, 0.0, -0.0]), SymOp([1.0 1.0 1.0; -1.0 0.0 0.0; 0.0 0.0 -1.0], [0.5, 0.0, -0.0]), SymOp([0.0 1.0 0.0; 0.0 0.0 1.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([0.0 0.0 -1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0, 0.0, -0.0]), SymOp([-1.0 -1.0 -1.0; 0.0 0.0 1.0; 1.0 0.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 -1.0 0.0; 1.0 1.0 1.0; 0.0 0.0 -1.0], [0.0, 0.5, -0.0]), SymOp([1.0 0.0 0.0; 0.0 1.0 0.0; -1.0 -1.0 -1.0], [0.0, 0.0, 0.5]), SymOp([1.0 1.0 1.0; 0.0 -1.0 0.0; -1.0 0.0 0.0], [0.5, 0.0, -0.0]), SymOp([0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.0, 0.0, -0.0]), SymOp([-1.0 0.0 0.0; 0.0 0.0 -1.0; 1.0 1.0 1.0], [0.0, 0.0, 0.5]), SymOp([0.0 0.0 1.0; -1.0 -1.0 -1.0; 0.0 1.0 0.0], [0.0, 0.5, -0.0])]
cell = Cell(lattice, base_positions, base_species, symops)
cbonds = canonical_bonds(cell, 2.)
dist2 = [distance(cell, b) for b=cbonds]

@assert dist1 ≈ dist2



### FCC lattice, primitive unit cell

lattice = Mat3(1, 1, 0,   1, 0, 1,   0, 1, 1) / 2
positions = [Vec3(0., 0, 0)]
species = [1]
cell = Cell(lattice, positions, species)

b1 = Bond{3}(1, 1, Vec3([1, 0, 0]))
b2 = Bond{3}(1, 1, Vec3([0, 1, 0]))
@assert S.is_equivalent_by_symmetry(cell, b1, b2)

cbonds = canonical_bonds(cell, 2.)
[distance(cell, b) for b=cbonds]

# Populate interactions for a random bond
bond = cbonds[4]
basis = S.basis_for_symmetry_allowed_couplings(cell, bond)
for x = basis
    display(clean_digits.(x, 4))
end
J = basis' * randn(length(basis))
S.verify_coupling_matrix(cell, bond, J)
bonds, Js = S.all_symmetry_related_interactions(cell, bond, J)

for (b, J) = zip(bonds, Js)
    display(b)
    display(J)
end

