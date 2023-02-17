@testitem "JET stability" begin
    using JET

    function test(mode)
        latvecs = lattice_vectors(1,1,2,90,90,90)
        crystal = Crystal(latvecs, [[0,0,0]])
        L = 2
        sys = System(crystal, (L,L,1), [SpinInfo(1, S=1)], mode)
        set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))
        polarize_spins!(sys, (0,0,1))

        # Test stability of LocalSampler
        sampler = LocalSampler(kT=0.2, propose=propose_flip)
        @test_opt step!(sys, sampler)

        # Test stability with mixed proposals
        propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
        sampler = LocalSampler(kT=0.2; propose)
        @test_opt step!(sys, sampler)

        # Test stability of Langevin
        langevin = Langevin(0.01, kT=0.2, λ=0.1)
        @test_opt step!(sys, langevin)

        # Test stability of ImplicitMidpoint. For some reason, isapprox() on
        # lists of SVectors is only type stable for Julia >= v1.9.
        if VERSION >= v"1.9-beta"
            integrator = ImplicitMidpoint(0.01)
            @test_opt step!(sys, integrator)
        end
    end

    test(:dipole)
    test(:SUN)
end
