@testitem "Type stability" begin
    using JET

    function test(mode)
        latvecs = lattice_vectors(1,1,2,90,90,90)
        crystal = Crystal(latvecs, [[0,0,0]])
        L = 2
        sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=2)], mode)

        # KB-TODO: Reenable these after optimizing expected_quadrupole and mul_quadrupole_matrices

        #=
        @test_opt energy(sys)
        
        sampler = LocalSampler(kT=0.2, propose=propose_flip)
        @test_opt step!(sys, sampler)

        propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
        sampler = LocalSampler(kT=0.2; propose)
        @test_opt step!(sys, sampler)

        langevin = Langevin(0.01, kT=0.2, λ=0.1)
        @test_opt step!(sys, langevin)

        integrator = ImplicitMidpoint(0.01)
        @test_opt step!(sys, integrator)
        =#
    end

    test(:dipole)
    test(:SUN)
end

# In the future this should be a static analysis, but for now we can explicitly
# run the functions. See https://github.com/aviatesk/JET.jl/issues/286
@testitem "Memory allocations" begin
    function test(mode)
        latvecs = lattice_vectors(1,1,2,90,90,90)
        crystal = Crystal(latvecs, [[0,0,0]])
        L = 2
        sys = System(crystal, (L,L,1), [SpinInfo(1, S=1, g=2)], mode)
        set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))
        polarize_spins!(sys, (0,0,1))

        # TODO: Diagonose possible @allocated bug. BenchmarkTools.@btime
        # suggests that this is actually zero-allocation.
        energy(sys)
        @test 16 >= @allocated energy(sys)

        propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
        sampler = LocalSampler(kT=0.2; propose)
        step!(sys, sampler)
        @test 0 == @allocated step!(sys, sampler)

        langevin = Langevin(0.01, kT=0.2, λ=0.1)
        step!(sys, langevin)
        @test 0 == @allocated step!(sys, langevin)

        integrator = ImplicitMidpoint(0.01)
        step!(sys, integrator)
        @test 0 == @allocated step!(sys, integrator)
    end

    test(:dipole)
    test(:SUN)
end
