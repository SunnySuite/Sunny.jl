@testitem "Type stability" begin
    using JET

    function test(mode)
        latvecs = lattice_vectors(1,1,2,90,90,90)
        crystal = Crystal(latvecs, [[0,0,0]])
        L = 2
        sys = System(crystal, [1 => Moment(s=1, g=2)], mode; dims=(L, L, 1))

        @test_opt energy(sys)
        
        sampler = LocalSampler(kT=0.2, propose=propose_flip)
        @test_opt step!(sys, sampler)

        propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
        sampler = LocalSampler(kT=0.2; propose)
        @test_opt step!(sys, sampler)

        langevin = Langevin(0.01; damping=0.1, kT=0.2)
        @test_opt step!(sys, langevin)

        integrator = ImplicitMidpoint(0.01)
        @test_opt step!(sys, integrator)
    end

    test(:dipole)
    test(:SUN)
end

# Cannot make this a static analysis (e.g. AllocCheck.jl) because there may be
# allocations when building error messages. The runtime checks below avoid these
# error paths.
@testitem "Memory allocations" begin
    function test(mode)
        latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
        crystal = Crystal(latvecs, [[0, 0, 0]])
        L = 2
        sys = System(crystal, [1 => Moment(s=1, g=2)], mode; dims=(L, L, 1))
        set_exchange!(sys, -1.0, Bond(1, 1, (1, 0, 0)))
        polarize_spins!(sys, [0, 0, 1])

        # Dynamic dispatch on System{N} allocates 16 bytes. Avoid this by
        # acquiring static type information.
        sys |> function(sys::System{N}) where N
            energy(sys)
            @test iszero(@allocated energy(sys))

            propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
            sampler = LocalSampler(kT=0.2; propose)
            step!(sys, sampler)
            @test iszero(@allocated step!(sys, sampler))

            langevin = Langevin(0.01; damping=0.1, kT=0.2)
            step!(sys, langevin)
            @test iszero(@allocated step!(sys, langevin))

            integrator = ImplicitMidpoint(0.01)
            step!(sys, integrator)
            @test iszero(@allocated step!(sys, integrator))
        end
    end

    test(:dipole)
    test(:SUN)
end
