@testitem "Reshaping entangled units" begin
    # Dimer oriented along x within a square cell, with weak inter-cell exchange.
    latvecs = [1.0 0 0; 0 1 0; 0 0 1]
    positions = [[0.0, 0, 0], [0.5, 0, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    # Bare system
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, 0.1, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, 0.1, Bond(2, 2, [0, 1, 0]))

    # Cannot reshape then entangle
    sys_sheared = resize_supercell(sys, (2, 1, 1))
    @test_throws "Entangle a single-cell system first, then reshape" Sunny.entangle_units(sys_sheared, [(1, 2)])

    # Entangle the original cell
    esys = Sunny.entangle_units(sys, [(1, 2)])
    set_q0(r) = for u in eachsite(r); set_coherent!(r, [0, 1/√2, -1/√2, 0], u); end
    set_q0(esys)
    E0 = energy_per_site(esys)

    # True when some unit has a member atom in a neighboring cell (a unit that
    # straddles the cell boundary of the reshaped system).
    straddles(r) = any(u -> any(id -> !iszero(id.cell_offset), u),
                       r.entanglement.contraction_info.inverse)

    # The `[-2 0 0; …]` shape flips the dimer's orientation so that a unit
    # straddles the cell boundary, exercising the `cell_offset` path in
    # `rebuild_entanglement!`/`member_site`. The others are ordinary reshapes.
    # Energy per site of a q=0 state is invariant under every commensurate reshape.
    for shape in ([-2 0 0; 0 1 0; 0 0 1],       # straddling
                  [1 1 0; -1 1 0; 0 0 1],       # shear
                  [3 0 0; 0 1 0; 0 0 1])        # resize ×3
        r = reshape_supercell(esys, shape)
        set_q0(r)
        @test energy_per_site(r) ≈ E0
    end
    @test straddles(reshape_supercell(esys, [-2 0 0; 0 1 0; 0 0 1]))

    # Full dynamics and coherent-state gradient through the straddling mapping.
    r = reshape_supercell(esys, [-2 0 0; 0 1 0; 0 0 1])
    set_field!(r, [0, 0, 10]); randomize_spins!(r); minimize_energy!(r)
    @test all(d -> isapprox(d[3], -1/2; atol=1e-6), r.entanglement.bare_system.dipoles)

    set_field!(r, [0.3, -0.7, 1.1]); randomize_spins!(r)
    Z0 = copy(r.coherents)
    ∇ = Sunny.energy_grad_coherents(r)
    maxerr = let ε = 1e-6, acc = 0.0
        for site in eachsite(r), c in 1:4, scale in (1.0 + 0im, im)
            δ = zeros(ComplexF64, 4); δ[c] = scale
            zn = copy(Z0); zn[site] = Z0[site] + ε*δ
            r.coherents .= zn; Sunny.set_expected_dipoles!(r)
            Ep = energy(r; check_normalization=false)
            zn[site] = Z0[site] - ε*δ
            r.coherents .= zn; Sunny.set_expected_dipoles!(r)
            Em = energy(r; check_normalization=false)
            acc = max(acc, abs((Ep - Em)/(2ε) - 2*real(δ' * ∇[site])))
        end
        acc
    end
    r.coherents .= Z0; Sunny.set_expected_dipoles!(r)
    @test maxerr < 1e-8

    # `repeat_periodically` agrees with the unreshaped energy per site.
    rep = repeat_periodically(esys, (2, 1, 1)); set_q0(rep)
    @test energy_per_site(rep) ≈ E0
end

# TODO: Add test with magnetic unit cell larger than a single unit (i.e. not q=0
# ordering).

@testitem "Dimer Tests" begin
    import LinearAlgebra: norm, I
    import Random

    J = 1.0
    J′ = 0.1
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0], [0.0, 0.5, 0.0]] 

    crystal = Crystal(latvecs, positions, 1; types = ["A", "B"])
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; seed=1)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, J′, Bond(2, 2, [1, 0, 0]))

    esys = Sunny.entangle_units(sys, [(1, 2)])
    bare = esys.entanglement.bare_system
    interactions = esys.interactions_union[1]

    # Test on-bond exchange
    onsite_operator = interactions.onsite
    S = spin_matrices(1/2)
    Sl, Su = to_product_space(S, S)
    onsite_ref = J * (Sl' * Su)
    @test onsite_operator ≈ onsite_ref

    # Test apparatus for setting coherent states from dipoles specification
    dipoles = [[0, 1/2, 0], [0, -1/2, 0]] # Dipoles specifying a dimer state
    cs = Sunny.coherent_state_from_dipoles(esys, dipoles, 1)
    set_coherent!(esys, cs, CartesianIndex(1,1,1,1))
    @test bare.dipoles[1,1,1,1][2] ≈ 1/2
    @test bare.dipoles[1,1,1,2][2] ≈ -1/2

    # Test external field. The Zeeman term is not folded into the onsite
    # operator; it lives in `extfield` (per unit) and couples to the unit's total
    # moment via the per-unit total-moment operator. The onsite operator is
    # unchanged.
    set_field!(esys, [0, 0, -10])
    @test esys.interactions_union[1].onsite ≈ onsite_ref
    @test esys.extfield[1] ≈ [0, 0, -10]
    @test bare.extfield[1, 1, 1, 1] ≈ [0, 0, -10]
    @test bare.extfield[1, 1, 1, 2] ≈ [0, 0, -10]

    # Test external field works in action
    randomize_spins!(esys)
    minimize_energy!(esys)
    @test bare.dipoles[1][3] ≈ 1/2
    @test bare.dipoles[2][3] ≈ 1/2
    set_field!(esys, [0, 0, 0])
    minimize_energy!(esys)
    @test norm(bare.dipoles[1]) < 1e-10
    @test norm(bare.dipoles[2]) < 1e-10

    # Test inter-bond exchange
    pc = Sunny.as_general_pair_coupling(interactions.pair[1], esys)
    Sl1, Sl2 = to_product_space(Sl, Sl)
    Su1, Su2 = to_product_space(Su, Su)
    bond_operator = zeros(ComplexF64, 16, 16)
    for (A, B) in pc.general.data
        bond_operator .+= kron(A, B)
    end
    bond_ref = J′*((Sl2' * Sl1) .+ (Su2' * Su1))
    @test bond_operator ≈ bond_ref

    # Test dispersion against analytical formula for antisymmetric channel.
    qs = [[0.2, 0.3, 0]]

    ω_ref(q, J, J′) = J*sqrt(1 + 2(J′/J) * cos(2π*q))
    ωs_analytical = ω_ref.([q[1] for q in qs], J, J′)

    set_field!(esys, [0, 0, 0])
    for unit in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], unit)
    end
    swt = SpinWaveTheory(esys; measure=Sunny.empty_measurespec(esys), regularization=0.0)
    disp = dispersion(swt, qs)
    ωs_numerical = disp[1,:]

    @test all(both -> isapprox(both[1], both[2]; atol=1e-12), zip(ωs_analytical, ωs_numerical))

    # Reshaped entangled system produces the correct intensities
    swt = SpinWaveTheory(esys; measure=ssf_perp(esys))
    res = intensities_bands(swt, qs)
    shape = [1 2 0; 0 1 0; 0 0 1]
    esys_sheared = reshape_supercell(esys, shape)
    swt_sheared = SpinWaveTheory(esys_sheared; measure=ssf_perp(esys_sheared))
    res_sheared = intensities_bands(swt_sheared, qs)
    @test res_sheared.disp ≈ res.disp atol=1e-11
    @test sum(res_sheared.data; dims=1) ≈ sum(res.data; dims=1) atol=1e-11

    # Test static structure factor is zero (dipolar sector)
    ssf = SampledCorrelationsStatic(esys; measure=ssf_trace(esys))
    add_sample!(ssf, esys)
    @test all(x -> isapprox(x, 0.0; atol=1e-12), ssf.parent.data)

    ### Golden test for classical dynamics ###

    esys = repeat_periodically(esys, (8, 1, 1))
    for site in eachsite(esys)
        set_coherent!(esys, [1, 0, 0, 0], site)
        i = site[1]
        set_field_at!(esys, [cos(i), sin(2i), cos(3i)], site)
    end

    energies = range(0, 2, 5)
    dt = 0.1
    measure = ssf_trace(esys)
    sc = SampledCorrelations(esys; dt, energies, measure)
    integrator = ImplicitMidpoint(dt)

    for _ in 1:100
        step!(esys, integrator)
    end
    add_sample!(sc, esys)
    res = intensities(sc, qs; energies, kT=0.05)

    @test res.data ≈ [0.6963636938867421; 4.833911994098928; 28.719089468355055; 52.61793172686838; 42.62520819431053;;]
end

@testitem "Entangled dipole-dipole" begin
    import LinearAlgebra: norm, I
    import Random

    units = Units(:meV, :angstrom)

    # Two spin-1/2 atoms per cell, well separated so the moments are physical
    # point dipoles. Intra-cell exchange makes them one entangled unit; a weak
    # inter-cell exchange gives a nontrivial ground state.
    latvecs = [3.0 0 0; 0 3 0; 0 0 6]
    positions = [[0, 0, 0], [0, 0.15, 0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    function make_entangled(; seed)
        sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; seed)
        set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]))
        set_exchange!(sys, 0.1, Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys, 0.1, Bond(2, 2, [1, 0, 0]))
        # `units` indexes the chemical cell; reshape the entangled system to size.
        esys = resize_supercell(Sunny.entangle_units(sys, [(1, 2)]), (2, 1, 1))
        enable_dipole_dipole!(esys, units.vacuum_permeability)
        return esys
    end

    # `enable_dipole_dipole!` builds Ewald on the physical bare system, not the
    # contracted system.
    esys = make_entangled(seed=1)
    @test isnothing(esys.ewald)
    @test !isnothing(esys.entanglement.bare_system.ewald)

    # The dipole-dipole energy (evaluated on the physical bare system) must
    # equal that of an equivalent ordinary SU(N) system with the same physical
    # spins, for a matching polarized state. Compare the Ewald term in
    # isolation, since `ord` carries no exchange couplings.
    ord = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; dims=(2, 1, 1), seed=1)
    enable_dipole_dipole!(ord, units.vacuum_permeability)
    polarize_spins!(ord, [0, 0, 1])
    for u in eachsite(esys)
        cs = Sunny.coherent_state_from_dipoles(esys, [[0, 0, 1/2], [0, 0, 1/2]], Sunny.to_atom(u))
        set_coherent!(esys, cs, u)
    end
    @test Sunny.ewald_energy(esys.entanglement.bare_system) ≈ Sunny.ewald_energy(ord)

    # The coherent-state gradient dE/dZ̄ (which includes the promoted Ewald term)
    # must match a finite difference.
    set_field!(esys, [0.3, -0.7, 1.1])
    randomize_spins!(esys)
    Z0 = copy(esys.coherents)
    ∇ = Sunny.energy_grad_coherents(esys)
    maxerr = let ε = 1e-6, acc = 0.0
        for site in eachsite(esys), c in 1:4, scale in (1.0 + 0im, im)
            δ = zeros(ComplexF64, 4); δ[c] = scale
            znew = copy(Z0); znew[site] = Z0[site] + ε*δ
            esys.coherents .= znew; Sunny.set_expected_dipoles!(esys)
            Ep = energy(esys; check_normalization=false)
            znew[site] = Z0[site] - ε*δ
            esys.coherents .= znew; Sunny.set_expected_dipoles!(esys)
            Em = energy(esys; check_normalization=false)
            acc = max(acc, abs((Ep - Em)/(2ε) - 2*real(δ' * ∇[site])))
        end
        acc
    end
    esys.coherents .= Z0; Sunny.set_expected_dipoles!(esys)
    @test maxerr < 1e-8

    # Ewald survives cloning and periodic repetition.
    @test !isnothing(Sunny.clone_system(esys).entanglement.bare_system.ewald)
    @test !isnothing(repeat_periodically(esys, (2, 1, 1)).entanglement.bare_system.ewald)

    # SpinWaveTheory supports entangled dipole-dipole: the real-space dipolar
    # coupling enters through the physical bare atoms (via the generalized Ewald
    # block over unit "parts"). A strong field polarizes the units into a stable
    # ground state; the dispersion is then finite, and it changes when the
    # dipole-dipole interaction is present.
    qs = [[0.1, 0, 0], [0.3, 0.2, 0]]
    set_field!(esys, [0, 0, 40])
    randomize_spins!(esys)
    minimize_energy!(esys)
    disp_dd = dispersion(SpinWaveTheory(esys; measure=nothing), qs)
    @test all(isfinite, disp_dd) && all(≥(-1e-6), disp_dd)

    let sys_nodd = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; seed=1)
        set_exchange!(sys_nodd, 1.0, Bond(1, 2, [0, 0, 0]))
        set_exchange!(sys_nodd, 0.1, Bond(1, 1, [1, 0, 0]))
        set_exchange!(sys_nodd, 0.1, Bond(2, 2, [1, 0, 0]))
        esys_nodd = resize_supercell(Sunny.entangle_units(sys_nodd, [(1, 2)]), (2, 1, 1))
        set_field!(esys_nodd, [0, 0, 40])
        randomize_spins!(esys_nodd)
        minimize_energy!(esys_nodd)
        disp_nodd = dispersion(SpinWaveTheory(esys_nodd; measure=nothing), qs)
        @test !isapprox(disp_dd, disp_nodd; atol=1e-6)
    end

    # LocalSampler still does not support long-range dipole-dipole for entangled units.
    let sampler = LocalSampler(kT=0.1, propose=Sunny.propose_flip)
        @test_throws "does not yet support long-range dipole-dipole" step!(esys, sampler)
    end
end

@testitem "Model parameters for entangled units" begin
    latvecs = [1 0 0; 0 1 0; 0 0 2]
    positions = [[0, 0, 0.0], [0, 0.5, 0.0]]
    crystal = Crystal(latvecs, positions, 1; types=["A", "B"])

    # A dimer (intra-unit :J) with an inter-unit exchange (:Jp) along x, so that
    # the entangled system's dispersion is pinned by both labeled parameters.
    sys = System(crystal, [1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)], :SUN; seed=1)
    set_exchange!(sys, 1.0, Bond(1, 2, [0, 0, 0]), :J => 1.0)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]), :Jp => 0.2)
    esys = Sunny.entangle_units(sys, [(1, 2)])
    for u in eachsite(esys)
        set_coherent!(esys, [0, 1/√2, -1/√2, 0], u)
    end

    qs = [[0.1, 0, 0], [0.25, 0, 0], [0.4, 0, 0]]
    dispersion_at(s) = let
        swt = SpinWaveTheory(s; measure=Sunny.empty_measurespec(s), regularization=0.0)
        dispersion(swt, qs)[1, :]
    end

    # `get_param` reads the labels off the physical (bare) system.
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]

    # `set_params!` mutably regenerates the contracted couplings, changing the
    # dispersion. A round-trip back to the original values reproduces it exactly.
    disp0 = dispersion_at(esys)
    set_params!(esys, [:J, :Jp], [1.5, 0.3])
    @test get_params(esys, [:J, :Jp]) == [1.5, 0.3]
    disp1 = dispersion_at(esys)
    @test !(disp1 ≈ disp0)
    set_params!(esys, [:J, :Jp], [1.0, 0.2])
    @test dispersion_at(esys) ≈ disp0

    # `make_loss_fn` operates on a clone: it sees the trial parameters, and the
    # original `esys` parameters are left untouched. The loss vanishes exactly at
    # the target parameters that generated `disp1`.
    loss = make_loss_fn(esys, [:J, :Jp]) do s
        return squared_error(disp1, dispersion_at(s))
    end
    @test loss([1.5, 0.3]) < 1e-20
    @test loss([1.0, 0.2]) > 1e-6
    @test get_params(esys, [:J, :Jp]) == [1.0, 0.2]
end

@testitem "Entangled Unit Intensity Scaling" begin
    latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
    positions = [[0, 0, 0], [0, 1/3, 0]]
    crystal = Crystal(latvecs, positions)

    J  = 1
    J′ = 0.12
    J2 = 0.1

    sys = System(crystal, [1 => Moment(s=1/2, g=2)], :SUN)
    set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, J′, Bond(1, 2, [0, -1, 0]));
    set_exchange!(sys, J2, Bond(1, 1, [1, 0, 0]));

    set_field!(sys, [0., 0., 10])

    esys = Sunny.entangle_units(sys, [(1, 2)])

    for s in (sys, esys)
        minimize_energy!(s)
        swt = SpinWaveTheory(s; measure=ssf_trace(s; apply_g=false))
        qs = q_space_path(crystal, [[0, 1, 0], [1/2, 1, 0], [1, 1, 0], [0, 0, 0]], 20)
        res = intensities_bands(swt, qs)
        @test sum(res.data[:,1]) ≈ 1
    end
end


# @testitem "Ba3Mn2O8 Dispersion and Golden Test" begin
# end
