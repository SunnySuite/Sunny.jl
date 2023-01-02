struct SFData 
    data        :: Array{ComplexF64, 7}                    # Raw SF data for 1st BZ 
    crystal     :: Crystal           
    Δω          :: Float64                                 # Energy step size
    idxinfo     :: SortedDict{CartesianIndex{2}, Int64}    # (α, β) to save from 𝒮^{αβ}(q, ω)
    site_infos  :: Vector{SiteInfo}                        # For form factor information
end

struct SFTrajectory{N}
    sys         :: SpinSystem{N}         # Clone system so original SpinSystem unaltered by trajectory calculation
    traj        :: Array{ComplexF64, 6}  # Trajectory buffer
    ops         :: Array{ComplexF64, 3}  # Operators corresponding to observables
    measperiod  :: Int                   # Steps to skip between saving observables (downsampling)
    gfactor     :: Bool
    dipolemode  :: Bool                  # Whether considering only dipoles 
    integrator  :: ImplicitMidpoint 
end

mutable struct StructureFactor{N}
    sfdata      :: SFData
    sftraj      :: SFTrajectory{N}
    nsamples    :: Int64
end


function StructureFactor(sys::SpinSystem;
    Δt = 0.1, numω = 100, ωmax = nothing, gfactor = true, ops = nothing, matrix_elems = nothing,
)
    sftraj = SFTrajectory(sys; Δt, numω, ωmax, ops, gfactor)
    sfdata = SFData(sys, sftraj; ops, matrix_elems)
    numsamps = 0

    return StructureFactor(sfdata, sftraj, numsamps)
end

function SFTrajectory(sys::SpinSystem{N}; 
    Δt = 0.1, numω = 100, ωmax = nothing, ops = nothing, gfactor = true,
) where N
    # Default to dipole expectation values if no observables have been given
    dipolemode = false 
    if isnothing(ops)
        dipolemode = true
        ops = zeros(ComplexF64, 0, 0, 3) # Placeholder with necessary information for consistent behavior later 
    else
        if N == 0 
            error("Structure Factor Error: Cannot provide matrices for observables for a dipolar `SpinSystem`")
        end
    end

    # Determine meas_period (downsampling factor)
    if isnothing(ωmax)
        measperiod = 1
    else
        @assert π/Δt > ωmax "Maximum ω with chosen step size is $(π/Δt). Choose smaller Δt or change ω_max."
        measperiod = floor(Int, π/(Δt * ωmax))
    end

    # Preallocation
    qa, qb, qc, ns = size(sys.dipoles)
    nops = size(ops, 3)
    traj = zeros(ComplexF64, nops, qa, qb, qc, ns, numω)
    integrator = ImplicitMidpoint(Δt)

    # Create a shallow copy of the spin system
    sys_new = SpinSystem(sys.crystal, sys.size, sys.hamiltonian,
        copy(sys.dipoles), copy(sys.coherents), sys.dipole_buffers, sys.coherent_buffers,
        sys.ℌ_buffer, sys.site_infos, sys.units, sys.rng)

    return SFTrajectory(sys_new, traj, ops, measperiod, gfactor, dipolemode, integrator)
end


function SFData(sys::SpinSystem, sftraj::SFTrajectory; 
    ops = nothing, matrix_elems = nothing,
)
    nops =  isnothing(ops) ? 3 : size(ops, 3) # Assume three observables (spin operators) if none are explicitly given
    numω = size(sftraj.traj, 6)

    # Save all matrix elements if subset isn't given
    if isnothing(matrix_elems)
        matrix_elems = []
        for i in 1:nops, j in i:nops
            push!(matrix_elems, (i, j))
        end
    end

    # Construct look-up table for matrix elements
    count = 1
    pairs = []
    for elem in matrix_elems
        α, β = elem
        α, β = α < β ? (α, β) : (β, α)  # Because SF is symmetric, only save diagonal and upper triangular
        push!(pairs, (α, β) => count)
        count += 1
    end
    pairs = map(i -> CartesianIndex(i.first) => i.second, pairs) # Convert to CartesianIndices
    idxinfo = SortedDict{CartesianIndex{2}, Int64}(pairs)

    qa, qb, qc, ns = size(sys.dipoles)
    data = zeros(ComplexF64, length(matrix_elems), qa, qb, qc, ns, ns, numω)
    Δω = 2π / (sftraj.integrator.Δt*sftraj.measperiod*numω)

    return SFData(data, sys.crystal, Δω, idxinfo, sys.site_infos) 
end


function Base.getindex(sfd::SFData, α, β, qa, qb, qc, l1, l2, ω)
    α, β = α < β ? (α, β) : (β, α)  # Because SF is symmetric, only save diagonal and upper triangular
    return sfd.data[sfd.idx_info[(α, β)], qa, qb, qc, l1, l2, ω]
end
Base.getindex(sf::StructureFactor, α, β, qa, qb, qc, l1, l2, ω) = sf.sfdata[α, β, qa, qb, qc, l1, l2, ω]


function calculate_structure_factor(sys::SpinSystem, sampler::LangevinSampler;
    ωmax=10.0, numω=100, numsamps=10, gfactor=true, Δt = nothing,
    ops = nothing, matrix_elems = nothing
)
    # Take a step size twice as large as the sampler step size if none explicitly given
    isnothing(Δt) && (Δt = 2sampler.integrator.Δt)

    sf = StructureFactor(sys; Δt, numω, ωmax, gfactor, ops, matrix_elems)
    for _ ∈ 1:numsamps
        sample!(sys, sampler)
        add_trajectory!(sf, sys)
    end

    return sf
end



include("SFUtils.jl")
include("Trajectories.jl")
include("FormFactor.jl")
include("ElementContraction.jl")
include("BasisReduction.jl")
include("Interpolation.jl")
include("PowderAveraging.jl")
include("DataRetrieval.jl")