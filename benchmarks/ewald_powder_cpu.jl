#!/usr/bin/env julia

# Run from the repository root as:
#   julia --project=. benchmarks/ewald_powder_cpu.jl [nq]

using LinearAlgebra
using Printf
using Random
using Sunny

function fibonacci_sphere(n)
    ϕ = (1 + √5) / 2
    return [begin
        z = 1 - 2(i - 0.5) / n
        r = √max(0, 1 - z^2)
        θ = 2π * (i / ϕ)
        Sunny.Vec3(r*cos(θ), r*sin(θ), z)
    end for i in 1:n]
end

function synthetic_dipolar_swt()
    latvecs = lattice_vectors(10.19, 10.19, 10.19, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]], 227)
    sys = System(cryst, [1 => Moment(s=7/2, g=2)], :dipole)
    sys = reshape_supercell(sys, primitive_cell(cryst))
    enable_dipole_dipole!(sys, Units(:K, :angstrom).vacuum_permeability)

    Random.seed!(1)
    randomize_spins!(sys)
    return SpinWaveTheory(sys; measure=nothing)
end

function benchmark_ewald_hamiltonian(; nq=10_000, radius=0.8)
    swt = synthetic_dipolar_swt()
    qs = radius .* fibonacci_sphere(nq)
    L = Sunny.nbands(swt)
    Hs = [zeros(ComplexF64, 2L, 2L) for _ in 1:Threads.maxthreadid()]
    checksums = zeros(Float64, Threads.maxthreadid())

    for q in qs[1:min(end, 32)]
        Sunny.swt_hamiltonian_dipole!(Hs[1], swt, q)
    end

    t = @elapsed Threads.@threads for iq in eachindex(qs)
        tid = Threads.threadid()
        H = Hs[tid]
        Sunny.swt_hamiltonian_dipole!(H, swt, qs[iq])
        checksums[tid] += real(tr(H))
    end

    @printf("synthetic_ewald_hamiltonian nq=%d threads=%d maxthreadid=%d\n",
            nq, Threads.nthreads(), Threads.maxthreadid())
    @printf("total %.6f per_q %.6f checksum %.12g\n", t, t/nq, sum(checksums))
end

nq = isempty(ARGS) ? 10_000 : parse(Int, ARGS[1])
benchmark_ewald_hamiltonian(; nq)
