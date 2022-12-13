struct SFData 
    data        :: Array{ComplexF64, 7}              # Raw SF data for 1st BZ (complex for off-diagonals)
    crystal     :: Crystal           
    Δω          :: Float64                           # Energy step size
    idx_info    :: SortedDict{Tuple{Int, Int}, Int}  # (α, β) to save from 𝒮^{αβ}(q, ω)
    site_infos  :: Vector{SiteInfo}                  # For form factor information
end

struct SFTrajectory
    sys         :: SpinSystem            # Clone system so original SpinSystem unaltered by trajectory calculation
    traj        :: Array{ComplexF64, 6}  # Trajectory buffer
    ops         :: Array{ComplexF64, 3}  # Operators corresponding to observables
    meas_period :: Int                   # Steps to skip between saving observables (downsampling)
    g_factor    :: Bool
    dipolemode  :: Bool                  # Whether considering only dipoles 
    integrator  :: ImplicitMidpoint 
end

mutable struct StructureFactor
    sfdata      :: SFData
    sftraj      :: SFTrajectory
    num_samples :: Int64
end

function clone_spin_system(sys::SpinSystem)
    (; 
        crystal, size, hamiltonian, dipoles, coherents, dipole_buffers, 
        coherent_buffers, ℌ_buffer, site_infos, consts, rng
    ) = sys
    dipoles_new = copy(dipoles)
    coherents_new = copy(coherents)
    return SpinSystem(crystal, size, hamiltonian, dipoles_new, coherents_new,
        dipole_buffers, coherent_buffers, ℌ_buffer, site_infos, consts, rng)
end


function StructureFactor(sys::SpinSystem{N};
    Δt = 0.1, num_ωs = 100, ω_max = nothing, g_factor = true,
    ops = nothing, matrix_elems = nothing
) where N

    sftraj = SFTrajectory(sys; Δt, num_ωs, ω_max, ops, g_factor)
    sfdata = SFData(sys, sftraj; ops, matrix_elems)

    return StructureFactor(sfdata, sftraj, 0)
end


function SFTrajectory(sys::SpinSystem{N}; 
    Δt = 0.1, num_ωs = 100, ω_max = nothing, ops = nothing, g_factor = true,
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
    if isnothing(ω_max)
        meas_period = 1
    else
        @assert π/Δt > ω_max "Maximum ω with chosen step size is $(π/Δt). Choose smaller Δt or change ω_max."
        meas_period = floor(Int, π/(Δt * ω_max))
    end

    # Preallocation
    qa, qb, qc, ns = size(sys.dipoles)
    nops = size(ops, 3)
    traj = zeros(ComplexF64, nops, qa, qb, qc, ns, num_ωs)
    integrator = ImplicitMidpoint(Δt)
    sys_new = clone_spin_system(sys)

    return SFTrajectory(sys_new, traj, ops, meas_period, g_factor, dipolemode, integrator)
end


function SFData(sys::SpinSystem, sftraj::SFTrajectory;
    ops = nothing, matrix_elems = nothing
)
    nops =  isnothing(ops) ? 3 : size(ops, 3) # Assume three observables (spin operators) if none are explicitly given
    num_ωs = size(sftraj.traj, 6)

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
    idx_info = SortedDict{Tuple{Int64, Int64}, Int64}(pairs)

    qa, qb, qc, ns = size(sys.dipoles)
    data = zeros(ComplexF64, length(matrix_elems), qa, qb, qc, ns, ns, num_ωs)
    Δω = 2π /(sftraj.integrator.Δt*sftraj.meas_period*num_ωs)

    return SFData(data, sys.crystal, Δω, idx_info, sys.site_infos) 
end


function Base.getindex(sfd::SFData, α, β, qa, qb, qc, l1, l2, ω)
    α, β = α < β ? (α, β) : (β, α)  # Because SF is symmetric, only save diagonal and upper triangular
    return sfd.data[sfd.idx_info[(α, β)], qa, qb, qc, l1, l2, ω]
end
Base.getindex(sf::StructureFactor, α, β, qa, qb, qc, l1, l2, ω) = sf.sfdata[α, β, qa, qb, qc, l1, l2, ω]


include("SFUtils.jl")
include("Trajectories.jl")
include("FormFactor.jl")
include("ElementContraction.jl")
include("BasisReduction.jl")
include("DataRetrieval.jl")
include("PowderAveraging.jl")