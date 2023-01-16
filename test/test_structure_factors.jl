@testitem "Sum Rule" begin

    function simple_model_sf(; mode)
        latsize = (4,4,4)
        J = 1.0
        cryst = Sunny.fcc_primitive_crystal()

        S = mode==:SUN ? 1/2 : 1
        κ = mode==:SUN ? 2 : 1
        sys = System(cryst, latsize, [SiteInfo(1; S)]; mode)
        sys.κs .= κ
        set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
        return sys
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
    
    for mode in (:SUN, :dipole)
        ops = if mode==:SUN
            N = 2
            S = Sunny.spin_matrices(N)
            ops = cat(S...; dims=3)
        else
            nothing
        end

        sys = simple_model_sf(; mode)
        sf = StructureFactor(sys; ωmax=10.0, gfactor=false, ops)
        thermalize_simple_model!(sys; kT=0.1)
        add_trajectory!(sf, sys)
        intensities = intensity_grid(sf; contraction=Trace(), negative_energies=true)

        @test isapprox(sum(intensities) / prod(sys.latsize), 1.0; atol=1e-12)
    end
end
