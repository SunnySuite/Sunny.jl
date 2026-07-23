################################################################################
# Model implementation
################################################################################

# Details about the Hamiltonian implemented below may be found in Bai et al.,
# "Hybridized quadrupolar excitations," Nature Physics 17
# (https://www.nature.com/articles/s41567-020-01110-1, https://arxiv.org/abs/2004.05623)


# Wave vectors for the three symmetry-equivalent ground states.
const q_gs = [
    [0, -1/4, 1/4],
    [1/4, 0, 1/4],
    [-1/4, 1/4, 1/4]
]

# Instantiate crystal in Sunny and restrict to magnetic (Fe) ions only.
function FeI2_crystal()
    a = b = 4.05012 
    c = 6.75214
    latvecs = lattice_vectors(a, b, c, 90, 90, 120)
    positions = [[0,0,0], [1/3, 2/3, 1/4], [2/3, 1/3, 3/4]]
    types = ["Fe", "I", "I"]
    FeI2 = Crystal(latvecs, positions; types)
    subcrystal(FeI2, "Fe")
end

# Generate a single unit cell in one of the ground states. Note that
# this performance a quick minimization process. (For speed, it would be
# better to save instantiated systems in each ground state.)
function FeI2_magnetic_unit_cell(; gs = 1, seed=nothing)
    gs_mats = [
        [1 0 0; 0 1 -2; 0 1 2],
        [1 0 2; 0 1 0; -1 0 2],
        [2 0 1; -1 1 0; -1 -1 1]
    ]
    cryst = FeI2_crystal()
    sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; seed, dims=(4,4,4))

    J1pm   = -0.236
    J1pmpm = -0.161
    J1zpm  = -0.261 
    J2pm   = 0.026
    J3pm   = 0.166
    J′0pm  = 0.037
    J′1pm  = 0.013
    J′2apm = 0.068
    J1zz   = -0.236
    J2zz   = 0.113
    J3zz   = 0.211
    J′0zz  = -0.036
    J′1zz  = 0.051
    J′2azz = 0.073
    J1xx = J1pm + J1pmpm 
    J1yy = J1pm - J1pmpm
    J1yz = J1zpm
    set_exchange!(sys, [J1xx 0.0 0.0; 0.0 J1yy J1yz; 0.0 J1yz J1zz], Bond(1,1,[1,0,0]))
    set_exchange!(sys, [J2pm 0.0 0.0; 0.0 J2pm 0.0; 0.0 0.0 J2zz], Bond(1,1,[1,2,0]))
    set_exchange!(sys, [J3pm 0.0 0.0; 0.0 J3pm 0.0; 0.0 0.0 J3zz], Bond(1,1,[2,0,0]))
    set_exchange!(sys, [J′0pm 0.0 0.0; 0.0 J′0pm 0.0; 0.0 0.0 J′0zz], Bond(1,1,[0,0,1]))
    set_exchange!(sys, [J′1pm 0.0 0.0; 0.0 J′1pm 0.0; 0.0 0.0 J′1zz], Bond(1,1,[1,0,1]))
    set_exchange!(sys, [J′2apm 0.0 0.0; 0.0 J′2apm 0.0; 0.0 0.0 J′2azz], Bond(1,1,[1,2,1]))
    D = 2.165
    S = spin_matrices(1)
    set_onsite_coupling!(sys, -D*S[3]^2, 1)

    sys = reshape_supercell(sys, gs_mats[gs])
    randomize_spins!(sys)
    minimize_energy!(sys)

    return sys, cryst
end

# Build a FeI2 system in one of the three ground states `gs`.
function FeI2_sys_and_cryst(dims=(4,4,4); gs=1, seed=nothing)
    sys, cryst = FeI2_magnetic_unit_cell(; gs, seed)
    resize_supercell(sys, dims), cryst
end