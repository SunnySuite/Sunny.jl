""" Functions for computing and manipulating structure factors """

## TODO: Structure factors should be OffsetArrays for a nice interface
## TODO: By default, return structure factors with basis indices still
##        present, then expose functions which do the phase-weighted
##        sum over the sublattices.

"""
    _fft_spin_traj(spin_traj, lattice; [plan, fft_space])

Takes in an array of spins (Vec3) of shape [B, D1, ..., Dd, T],
 with D1 ... Dd being the spatial dimensions, B the sublattice index,
 and T the time index.
Computes and returns an array of the shape [3, D1, ..., Dd, T],
 holding spatial and temporal fourier transforms S^α(q, ω), and the 
 phase-weighted sum already performed over the basis index. To get the
 contribution to the structure factor matrix, you want to then perform:
    𝒮^αβ(q, ω) = S^α(q, ω) * S^β(q, ω)∗
"""
function _fft_spin_traj(
    spin_traj::Array{Vec3}, lattice::Lattice;
    plan::Union{Nothing, FFTW.cFFTWPlan}=nothing,
    fft_space::Union{Nothing, Array{ComplexF64}}=nothing,
)
    @assert size(spin_traj)[1:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"
    if !isnothing(fft_space)
        @assert size(fft_space) == tuplejoin(3, size(spin_traj)) "fft_space size not compatible with spin_traj size"
    end

    # Reinterpret array to add the spin dimension explicitly
    # Now of shape [3, B, D1, ..., Dd, T]
    spin_traj = _reinterpret_from_spin_array(spin_traj)

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        fft_space = fft(spin_traj, 3:ndims(spin_traj))
    elseif isnothing(fft_space)
        fft_space = plan * spin_traj
    else
        mul!(fft_space, plan, spin_traj)
    end

    # Multiply each sublattice by the appropriate phase factors
    #    S_b(q, ω) → exp(-i b ⋅ q) S_b(q, ω)
    wave_vectors = Iterators.product((fftfreq(size(spin_traj, 2+d)) for d in 1:D)...)
    for (q_idx, q) in zip(eachcellindex(lattice), wave_vectors)
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            fft_space[1:end, b_idx, q_idx, 1:end] .*= exp(-im * (b ⋅ q))
        end
    end

    # Sum over sublattices b (now with phase factors)
    fft_spins = dropdims(sum(fft_space; dims=2); dims=2)

    return fft_spins
end

"""
    _fft_spin_traj!(fft_spins, spin_traj, lattice; [plan])

Doubly in-place version of `_fft_spin_traj`. This form requires the spin
 trajectory to be in ComplexF64 eltype to allow in-place FFTs.
`spin_traj` will be updated with its FFT, and `fft_spins` will be updated
 with the appropriate phase-weighted sum across basis sites.
"""
function _fft_spin_traj!(
    fft_spins::Array{ComplexF64}, spin_traj::Array{ComplexF64}, lattice::Lattice{D};
    plan::Union{Nothing, FFTW.cFFTWPlan}=nothing,
) where {D}
    @assert size(spin_traj)[2:end-1] == size(lattice) "Spin trajectory array size not compatible with lattice"
    @assert size(fft_spins) == tuplejoin(3, size(spin_traj)[3:end]) "fft_spins size not compatible with spin_traj size"

    recip = gen_reciprocal(lattice)

    # FFT along the D spatial indices, and the T time index
    if isnothing(plan)
        fft!(spin_traj, 3:ndims(spin_traj))
    else
        spin_traj = plan * spin_traj
    end

    _phase_weight_basis!(fft_spins, spin_traj, lattice)
    return fft_spins
end

"""Weight each FFT amplitude by the basis-dependent phase, sum over basis sites,
    and write the results to `result`. Sizes should be of the form:
   `result`:       3 × L1 × ⋯ × LN × T
   `spins_ft`: 3 × B × L1 × ⋯ × LN × T
"""
function _phase_weight_basis!(result::Array{ComplexF64}, spin_traj_ft::Array{ComplexF64}, lattice::Lattice{D}) where {D}
    @assert size(result, 1) == 3 "First dimension of `result` should be length 3"
    @assert size(spin_traj_ft, 1) == 3 "First dimension of `spin_traj_ft` should be length 3"
    @assert size(result)[2:end] == size(spin_traj_ft)[3:end] "Sizes of `result` and `spin_traj_ft` don't agree!"
    recip = gen_reciprocal(lattice)
    # Total number of timesteps / frequency bins
    T = size(spin_traj_ft, ndims(spin_traj_ft))

    # Multiply each sublattice by the appropriate phase factors
    #    S_b(q, ω) → exp(-i b ⋅ q) S_b(q, ω)
    for q_idx in eachcellindex(lattice)
        q = recip.lat_vecs * Vec3((Tuple(q_idx) .- 1) ./ lattice.size)
        for (b_idx, b) in enumerate(lattice.basis_vecs)
            # WARNING: Cannot replace T with `end` due to Julia's
            #  mishandling of CartesianIndex and `end`.
            spin_traj_ft[1:end, b_idx, q_idx, 1:T] .*= exp(-im * (b ⋅ q))
        end
    end

    # Sum over sublattices b (now with phase factors)
    # Need to add a singleton dimension to struct_factor matching
    #  up with the basis index of fft_spins
    result = reshape(result, 3, 1, size(result)[2:end]...)
    sum!(result, spin_traj_ft)

    return dropdims(result; dims=2)
end

""" Applies the dipole form factor, reducing the structure factor tensor to the
     observable quantities. I.e. performs the contraction:
        ∑_αβ (δ_αβ - 𝐪̂_α 𝐪̂_β) 𝒮^αβ(𝐪, ω)
    `struct_factor` should be of size 3 × 3 × L1 × ⋯ × LN × T
    Returns a real array of size              L1 × ⋯ × LN × T
"""
function dipole_form_factor(struct_factor::Array{ComplexF64}, lattice::Lattice{D}) where {D}
    recip = gen_reciprocal(lattice)
    T = size(struct_factor, ndims(struct_factor))
    result = zeros(Float64, size(struct_factor)[3:end])
    for q_idx in eachcellindex(lattice)
        q = recip.lat_vecs * Vec3((Tuple(q_idx) .- 1) ./ lattice.size)
        q = q / norm(q)
        dip_factor = reshape(I(D) - q * q', 3, 3, 1)
        for t in 1:T
            result[q_idx, t] = real(dot(dip_factor, struct_factor[:, :, q_idx, t]))
        end
    end
    return result
end

function _plan_spintraj_fft(spin_traj::Array{Vec3})
    spin_traj = _reinterpret_from_spin_array(spin_traj)
    return plan_fft(spin_traj, 3:ndims(spin_traj))
end

function _plan_spintraj_fft!(spin_traj::Array{ComplexF64})
    return plan_fft!(spin_traj, 3:ndims(spin_traj))
end

"""
    full_structure_factor(sys, kt; nsamples, langevinΔt, langevin_steps, measureΔt,
                          measure_steps, collect_steps, α, verbose)

Measures the full structure factor tensor of a spin system.
Returns 𝒮^αβ = ⟨S^α(q, ω) S^β(q, ω)∗⟩, which is an array of shape
    [3, 3, D1, ..., Dd, T]
"""
function full_structure_factor(
    sys::SpinSystem{D}, kT::Float64; nsamples::Int=10, langevinΔt::Float64=0.01,
    langevin_steps::Int=10000, measureΔt::Float64=0.01, measure_steps::Int=100,
    collect_steps::Int=10, α::Float64=0.1, verbose::Bool=false,
) where {D}
    num_snaps = 1 + div(measure_steps - 1, collect_steps)
    nb = nbasis(sys.lattice)
    spin_traj = zeros(ComplexF64, 3, size(sys)..., num_snaps)
    fft_spins = zeros(ComplexF64, 3, size(sys)[2:end]..., num_snaps)
    struct_factor = zeros(ComplexF64, 3, 3, size(sys)[2:end]..., num_snaps)
    plan = _plan_spintraj_fft!(spin_traj)
    integrator = HeunP(sys)
    integratorL = LangevinHeunP(sys, kT, α)

    println("Beginning thermalization...")

    # Thermalize the system for five times langevin_steps
    @showprogress 1 "Initial thermalization " for _ in 1:5*langevin_steps
        evolve!(integratorL, langevinΔt)
    end

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    for n in 1:nsamples
        printstyled("Obtaining sample $n\n", bold=true)

        # Evolve under Langevin dynamics to sample a news state
        @showprogress 1 "Langevin sampling new state " for _ in 1:langevin_steps
            evolve!(integratorL, langevinΔt)
        end

        # Evolve at fixed energy to collect dynamics info
        selectdim(spin_traj, ndims(spin_traj), 1) .= _reinterpret_from_spin_array(sys.sites)
        @showprogress 1 "Measuring structure factor " for nsnap in 2:num_snaps
            for _ in 1:collect_steps
                evolve!(integrator, measureΔt)
            end
            selectdim(spin_traj, ndims(spin_traj), nsnap) .= _reinterpret_from_spin_array(sys.sites)
        end

        _fft_spin_traj!(fft_spins, spin_traj, sys.lattice; plan=plan)
        FSα = reshape(fft_spins, 1, size(fft_spins)...)
        FSβ = reshape(fft_spins, size(fft_spins, 1), 1, size(fft_spins)[2:end]...)
        @. struct_factor += FSα * conj(FSβ)
    end

    struct_factor ./= nsamples

    return struct_factor
end

"""
    diag_structure_factor(sys, kt; nsamples, langevinΔt, langevin_steps, measureΔt,
                          measure_steps, collect_steps, α, verbose)

Measures the diagonal elements of the structure factor of a spin system.
Returns 𝒮^α = ⟨S^α(q, ω) S^α(q, ω)∗⟩, which is an array of shape
    [3, D1, ..., Dd, T]
"""
function diag_structure_factor(
    sys::SpinSystem{D}, kT::Float64; nsamples::Int=10, langevinΔt::Float64=0.01,
    langevin_steps::Int=10000, measureΔt::Float64=0.01, measure_steps::Int=100,
    collect_steps::Int=10, α::Float64=0.1, verbose::Bool=false,
) where {D}
    num_snaps = 1 + div(measure_steps - 1, collect_steps)
    spin_traj = zeros(ComplexF64, 3, size(sys)..., num_snaps)
    fft_spins = zeros(ComplexF64, 3, size(sys)[2:end]..., num_snaps)
    struct_factor = zeros(Float64, 3, size(sys)[2:end]..., num_snaps)
    plan = _plan_spintraj_fft!(spin_traj)
    integrator = HeunP(sys)
    integratorL = LangevinHeunP(sys, kT, α)

    # Thermalize the system for ten times langevin_steps
    @showprogress 1 "Initial thermalization " for _ in 1:10*langevin_steps
        evolve!(integratorL, langevinΔt)
    end

    if verbose
        println("Done thermalizing. Beginning measurements...")
    end

    for n in 1:nsamples
        printstyled("Obtaining sample $n\n", bold=true)

        # Evolve under Langevin dynamics to sample a new state
        @showprogress 1 "Langevin sampling new state " for _ in 1:langevin_steps
            evolve!(integratorL, langevinΔt)
        end

        # Evolve at fixed energy to collect dynamics info
        selectdim(spin_traj, ndims(spin_traj), 1) .= _reinterpret_from_spin_array(sys.sites)
        @showprogress 1 "Measuring new state structure factor " for nsnap in 2:num_snaps
            for _ in 1:collect_steps
                evolve!(integrator, measureΔt)
            end
            selectdim(spin_traj, ndims(spin_traj), nsnap) .= _reinterpret_from_spin_array(sys.sites)
        end

        _fft_spin_traj!(fft_spins, spin_traj, sys.lattice; plan=plan)
        @. struct_factor += real(fft_spins * conj(fft_spins))
    end

    struct_factor ./= nsamples

    return struct_factor
end