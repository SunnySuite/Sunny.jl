using Sunny
import Sunny: Vec3
using GLMakie

function path_points(points; points_per_leg = 10)
    legs = []
    for i ∈ 1:length(points)-1
        leg = []
        p1, p2 = points[i], points[i+1]
        for n in 1:points_per_leg
            push!(leg, Vec3((1 - (n-1)/points_per_leg)*p1 + (n-1)*p2/points_per_leg))
        end
        push!(legs, leg)
    end
    push!(legs[end], Vec3(points[end]))
    return vcat(legs...)
end

function path(sf, points; points_per_leg = nothing)
    points_per_leg = isnothing(points_per_leg) ? round(Int, size(sf.sfdata.data)[2]/2) : points_per_leg
    println("points_per_leg", points_per_leg)
    path = path_points(points; points_per_leg)
    nω = size(sf.sfdata.data, 7)
    hnω = round(Int, nω/2)
    ωs = (0:size(sf.sfdata.data)[end]-hnω) .* sf.sfdata.Δω
    intensities = zeros(length(path), length(ωs))
    for (j, ω) in enumerate(ωs)
        for (i, q) in enumerate(path)
            intensities[i, j] = get_intensity(sf, q, ω;
                c2q_temp=2.0, contraction = trace)
        end
    end
    return intensities, ωs
end


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
    N = 2
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
    ω_max = 5.5        # Maximum frequency to resolve. Sunny will downsample
    num_ωs = 200       # Total number of frequencies we'd like to resolve

    sf = StructureFactor(sys; Δt, num_ωs, ω_max)
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

@time begin
    nsamples = 3
    for _ in 1:nsamples
        sample!(sys, sampler)
        add_trajectory!(sf, sys)
    end
end

begin
    ω = 3.0
    q = [π, π, 0.0]
    @time intensity = get_intensity(sf, q)
end

# Extract path
begin
    qz = 0.0 
    points = [[0.0, 0.0, qz], [π, 0.0, qz], [π, π, qz], [0.0, 0.0, qz]]
    intensities, ωs = path(sf, points; points_per_leg=100)
end;

begin
    heatmap(1:size(intensities, 1), ωs, intensities; colorrange=(0.0, 1.0))
end


begin
    @time intensities_all = Sunny.get_intensities(sf)
    size(intensities_all)
end


begin
    qs = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1], [π, π, π]]
    # qs = Sunny.nearest_q.(Ref(sf.sfdata), qs)
    # qis = [q[2] for q in qs]
    # sort(qis)
end

begin
    @time an_intensity = Sunny.get_intensities(sf, qs[1])
    # @time some_intensities = Sunny.get_intensities(sf, qs)
end;