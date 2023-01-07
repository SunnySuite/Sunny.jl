@testitem "Sum Rule" begin

    function simple_model_sf(; SUN)
        dims = (4,4,4)
        J = 1.0
        cryst = Sunny.fcc_primitive_crystal()
        interactions = [heisenberg(J, Bond(1, 1, [1, 0, 0]))]
        if SUN
            sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; S=1/2)]; SUN)
            sys.κs .= 2
            return sys
        else
            return SpinSystem(cryst, interactions, dims, [SiteInfo(1; S=1)]; SUN)
        end
    end

    function thermalize_simple_model!(sys; kT)
        Δt = 0.05  # Time step for thermalization
        λ  = 0.1
        nsteps = 1000  # Number of steps between MC samples
        integrator = LangevinHeunP(kT, λ, Δt)
        for _ ∈ 1:nsteps
            step!(sys, integrator)
        end
    end
    
    for SUN in (true, false)
        ops = if SUN
            N = 2
            S = Sunny.spin_matrices(N)
            ops = cat(S...; dims=3)
        else
            nothing
        end

        sys = simple_model_sf(; SUN)
        sf = StructureFactor(sys; ωmax=10.0, gfactor=false, ops)
        thermalize_simple_model!(sys; kT=0.1)
        add_trajectory!(sf, sys)
        intensities = intensity_grid(sf; contraction=Trace(), negative_energies=true)

        @test isapprox(sum(intensities) / prod(sys.size), 1.0; atol=1e-12)
    end
end
