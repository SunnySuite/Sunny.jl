@testitem "Sum Rule" begin

    function simple_model_sf(; dims = (4, 4, 4), N=0)
        J = 1.0
        cryst = Sunny.fcc_primitive_crystal()
        interactions = [heisenberg(J, Bond(1, 1, [1, 0, 0]))]
        spin_rescaling = N == 0 ? 1.0 : 2.0/(N-1)
        sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N, spin_rescaling)])
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
    
    # Trial parameters
    ωmax = 10.0
    dims = (4, 4, 4)
    gfactor = false
    Ns = [0, 2]

    for N in Ns
        ops = if N == 2
            Ss = Sunny.spin_matrices(2)
            ops = zeros(ComplexF64, 2, 2, 3)
            for i ∈ 1:3
                ops[:,:,i] = Ss[i]
            end 
            ops
        else
            nothing
        end

        sys = simple_model_sf(; N)
        sf = StructureFactor(sys; ωmax, gfactor, ops)
        thermalize_simple_model!(sys; kT=0.1)
        add_trajectory!(sf, sys)
        intensities = intensity_grid(sf; contraction=Trace(), negative_energies=true)

        @test isapprox(sum(intensities) / prod(dims), 1.0; atol=1e-12)
    end
end