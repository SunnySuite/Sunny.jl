@testset "Structure Factors" begin

function diamond_heisenberg_model(; 
    ff_elem=nothing, ff_lande=nothing, spin_rescaling, J
)
    crystal = Sunny.diamond_crystal()
    interactions = [
        heisenberg(J, Bond(1, 3, [0,0,0])),
    ]
    dims = (8, 8, 8)
    site_infos = [SiteInfo(1; spin_rescaling, ff_elem, ff_lande)]
    sys = SpinSystem(crystal, interactions, dims, site_infos) 
    rand!(sys)

    return sys
end

#= There are many different paths through the structure factor code.
This test simply tries them all out to make sure nothing is obviously
broken. It does not text for the correctness of the results, just for the
absence of errors. =#
function test_structure_factors_are_operational()
    dipole_factor = [true, false]
    reduce_basis = [true, false]
    ff_elem = [nothing, "Fe2"]
    ff_lande = [nothing, 3/2]

    for (dipole_factor, reduce_basis, ff_elem, ff_lande) in Base.Iterators.product(
         dipole_factor, reduce_basis, ff_elem, ff_lande
    )
        (!isnothing(ff_lande) && isnothing(ff_elem)) && continue

        J = Sunny.meV_per_K * 7.5413       
        spin_rescaling = 3/2
        sys = diamond_heisenberg_model(; ff_elem, ff_lande, spin_rescaling, J)

        Δt = 0.02 / (spin_rescaling^2 * J) 
        kT = Sunny.meV_per_K * 2. 
        α  = 0.1
        nsteps = 1  
        sampler = LangevinSampler(sys, kT, α, Δt, nsteps)

        dynΔt = Δt * 10
        omega_max = 5.5 
        num_omegas = 2

        dynamic_structure_factor(
            sys, sampler; nsamples=2, dt = dynΔt, omega_max, num_omegas,
            bz_size=(1,1,2), thermalize=1, verbose=false,
            reduce_basis, dipole_factor,
        )

        # Did we make it?
        @test true
    end

end

## Running these tests takes about 30 seconds. Not tested by default.
## Uncomment to run the test.
test_structure_factors_are_operational()

# TODO: Add test for correctness of calculation

end