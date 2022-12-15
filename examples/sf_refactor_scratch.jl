using Sunny
import Sunny: Vec3
using LinearAlgebra: norm
using GLMakie



function diamond_model(; dims = (2,2,2), SUN=true)
    crystal = Sunny.diamond_crystal()
    J = Sunny.CONSTS_meV.kB * 7.5413        # Units of meV
    interactions = [
        heisenberg(J, Bond(1, 3, [0,0,0])),
    ]
    N = SUN ? 2 : 0
    spin_rescaling = N == 2 ? 3 : 3/2
    sys = SpinSystem(crystal, interactions, dims, [SiteInfo(1; spin_rescaling, N)])
    rand!(sys)
    sys
end

################################################################################
# Test trajectory functions
################################################################################
begin
    sys = diamond_model()
    Δt = 0.05
    integrator = ImplicitMidpoint(Δt)
    Ss = Sunny.spin_matrices(2)
    ops = zeros(ComplexF64, 2, 2, 3)
    for i ∈ 1:3
        ops[:,:,i] = Ss[i]
    end 
    num_snaps = 1000

    begin
        # expectation_traj = Sunny.expectation_trajectory(integrator, num_snaps, Δt, ops)
        @time expectation_traj = Sunny.expectation_trajectory(sys, integrator, num_snaps, ops)
    end

    begin
        # dip_traj = Sunny.dipole_trajectory(integrator, num_snaps, Δt)
        @time dip_traj = Sunny.dipole_trajectory(sys, integrator, num_snaps; g_factor=false)
    end
end;

begin
    fig = lines(1:size(dip_traj, 6), real.(expectation_traj[1,1,1,1,1,:]))
    lines!((size(dip_traj, 6)+1):(2*size(dip_traj, 6)), real.(dip_traj[1,1,1,1,1,:]))
    fig
end




################################################################################
# Test StructureFactor and add_sample! 
################################################################################
begin
    sys = diamond_model(; SUN=true)
    Ss = Sunny.spin_matrices(2)
    ops = zeros(ComplexF64, 2, 2, 3)
    for i ∈ 1:3
        ops[:,:,i] = Ss[i]
    end 
    Δt = 0.05

    sftraj = Sunny.SFTrajectory(sys;
        ops, Δt, ω_max = 10.0, g_factor = false
    ) 
    sfdata = Sunny.SFData(sys, sftraj)

    @time new_trajectory!(sftraj, sys)
    @time accum_trajectory!(sfdata, sftraj, 1)
end;

begin
    sf = StructureFactor(sfdata, sftraj, 1)
end;

begin
    @time add_trajectory!(sf, sys)
end

val = Sunny.get_intensity(sf.sfdata, Sunny.Vec3(π, π, π), 0.0; c2q_temp=0.01)



################################################################################
# Test full calculation with basic tools available
################################################################################
begin
    dims = (8, 8, 8)
    N = 0
    SUN = N == 2 ? true : false
    spin_rescaling = N == 2 ? 3 : 3/2
    sys = diamond_model(; dims, SUN)

    spin_rescaling = 3/2
    J = Sunny.CONSTS_meV.kB * 7.5413        # Units of meV
    Δt_therm = 0.05 / (spin_rescaling^2 * J)     # Units of 1/meV
    kT = Sunny.CONSTS_meV.kB * 2. # Units of meV
    λ  = 0.1
    nsteps = 1000  # Number of steps between MC samples
    integrator = LangevinHeunP(kT, λ, Δt_therm)
    sampler = LangevinSampler(integrator, nsteps)

    Δt = 2*Δt_therm
    ωmax = 5.5        # Maximum frequency to resolve. Sunny will downsample
    numω = 200       # Total number of frequencies we'd like to resolve
    gfactor = false

    sf = StructureFactor(sys; Δt, numω, ωmax, gfactor)
end;

# Thermalize
@time begin
    for _ ∈ 1:10    
        sample!(sys, sampler)
    end
end


# Check near ground state
begin
    plot_spins(sys)
end

# All in one
@time begin
    nsamples = 3
    for _ in 1:nsamples
        sample!(sys, sampler)
        add_trajectory!(sf, sys)
    end
end

begin
    q = [π, π, 0.0]
    @time intensities = get_intensity(sf, q; interp=Sunny.NoInterp(), contraction=trace, c2q_temp=kT)
    @time intensity = Sunny.get_static_intensity(sf, q; interp=Sunny.NoInterp(), contraction=trace, c2q_temp=kT)
end;

# Extract path
begin
    qz = 0.0 
    points = [[0.0, 0.0, qz], [π, 0.0, qz], [π, π, qz], [0.0, 0.0, qz]]
    ωs = ωvals(sf)
    @time intensities = path(sf, points; 
        interp=Sunny.NoInterp(), contraction=depolarize, density=10, temp=kT,
    )
    heatmap(1:size(intensities, 1), ωs, intensities; colorrange=(0.0, 1.00))
end

@time begin
    intensity_grid = Sunny.get_intensity_grid(sf; contraction=Sunny.trace)
    sum(intensity_grid) / 8^4
end

begin
    qs = [(0.0, 0.0, 0.0) (0.0, 0.0, 0.1);
          (0.0, π, 0.0)   (π, π, π)]
    # qs = Sunny.nearest_q.(Ref(sf.sfdata), qs)
    # qis = [q[2] for q in qs]
    # sort(qis)
end

begin
    # @time an_intensity = Sunny.get_intensity(sf, qs[1])
    @time some_intensities = Sunny.get_intensities(sf, qs)
    @time static_intensities = Sunny.get_static_intensities(sf, qs)
end

begin
    (; intensities, qpoints) = Sunny.static_slice(sf, [0, 0, 0.0], [2π, 2π, 0.0], 0; 
        density=5, index_labels=true
    )
    heatmap(intensities)
end



# Check sum rule on simpler model
begin
    J = 4.0
    cryst = Sunny.fcc_primitive_crystal()
    interactions = [heisenberg(J, Bond(1, 1, [1, 0, 0]))]
    dims = (20, 20, 4)
    sys = SpinSystem(cryst, interactions, dims, [SiteInfo(1; N=0, spin_rescaling=1.0)])

    Δt_therm = 0.05/J
    kT = Sunny.CONSTS_meV.kB * 1. # Units of meV
    λ  = 0.1
    nsteps = 1000  # Number of steps between MC samples
    integrator = LangevinHeunP(kT, λ, Δt_therm)
    sampler = LangevinSampler(integrator, nsteps)

    Δt = 2*Δt_therm
    ωmax = 5.5        # Maximum frequency to resolve. Sunny will downsample
    numωs = 200       # Total number of frequencies we'd like to resolve
    g_factor = false

    sf = StructureFactor(sys; Δt, num_ωs, ω_max, g_factor)
end;




# Thermalize
@time begin
    for _ ∈ 1:25    
        sample!(sys, sampler)
    end
end

# Check near ground state
begin
    plot_spins(sys)
end

@time begin
    nsamples = 1
    for _ in 1:nsamples
        sample!(sys, sampler)
        add_trajectory!(sf, sys)
    end
end


begin
    qz = 0.0 
    points = [[0, 0, 0], [2π, 2π, 0]]
    intensities, ωs = path(sf, points; 
        interp=Sunny.NoInterp(), contraction=trace, density=10, c2q_temp=nothing
    )
    heatmap(1:size(intensities, 1), ωs, intensities; colorrange=(0.0, 0.01))
end

begin
    intensity_grid = Sunny.get_intensity_grid(sf; contraction=Sunny.trace)
    sum(intensity_grid) / (20*20*4)
end


################################################################################
# Try all in one version
################################################################################
begin
    dims = (8, 8, 8)
    N = 0
    SUN = N == 2 ? true : false
    spin_rescaling = N == 2 ? 3 : 3/2
    sys = diamond_model(; dims, SUN)

    spin_rescaling = 3/2
    J = Sunny.CONSTS_meV.kB * 7.5413        # Units of meV
    Δt_therm = 0.05 / (spin_rescaling^2 * J)     # Units of 1/meV
    kT = Sunny.CONSTS_meV.kB * 2. # Units of meV
    λ  = 0.1
    nsteps = 1000  # Number of steps between MC samples
    integrator = LangevinHeunP(kT, λ, Δt_therm)
    sampler = LangevinSampler(integrator, nsteps)

    ωmax = 5.5        # Maximum frequency to resolve. Sunny will downsample
    numω = 200       # Total number of frequencies we'd like to resolve
    gfactor = true
end;

# Thermalize
@time begin
    for _ ∈ 1:10    
        sample!(sys, sampler)
    end
end


@time begin
    sf = calculate_structure_factor(sys, sampler; numsamps=3, ωmax, numω)
end;

begin
    qz = 0.0 
    points = [[0.0, 0.0, qz], [π, 0.0, qz], [π, π, qz], [0.0, 0.0, qz]]
    ωs = ωvals(sf)
    hω = div(length(ωs), 2) + 1
    ωs = ωs[1:hω]
    intensities = path(sf, points; 
        interp=Sunny.NoInterp(), contraction=Sunny.Trace(), density=50, temp=kT
    )
    heatmap(1:size(intensities, 1), ωs, intensities; colorrange=(0.0, 1.0))
end

begin
    intensity_grid = Sunny.intensity_grid(sf; contraction=Sunny.Element((1,1)))
    sum(intensity_grid) / 8^4
end