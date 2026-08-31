# Standalone benchmark for building the q-dependent Ewald matrix A[cell, i, j].
#
# `precompute_dipole_ewald_at_wavevector` is planned by default: each call builds
# a `DipoleEwaldPlan` (hoisting all q-independent work) and materializes A. When
# sweeping many wavevectors, the plan should be built once and reused, which this
# benchmark quantifies.
#
# Run from the Sunny package environment:
#   include("benchmarks/ewald_qplan.jl")

using Sunny, LinearAlgebra, BenchmarkTools

# A crystal with several atoms per cell exercises the na² pair loop.
function test_crystal(; na=8)
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    positions = [[i/na, i/(2na), i/(3na)] for i in 0:na-1]
    Crystal(latvecs, positions, 1)
end

function run(; na=8, dims=(2,2,2), nq=200)
    cryst = test_crystal(; na)
    demag = Sunny.Mat3(I)  # trace-1 demag tensor for vacuum background
    qs = [Sunny.Vec3(0.1rand(), 0.1rand(), 0.1rand()) for _ in 1:nq]

    plan = Sunny.DipoleEwaldPlan(cryst, dims, demag)
    println("na=$na  dims=$dims  nq=$nq")
    println("plan build: ", @belapsed(Sunny.DipoleEwaldPlan($cryst, $dims, $demag)), " s")

    println("\n-- rebuild plan every q (public entry point) --")
    @btime for q in $qs
        Sunny.precompute_dipole_ewald_at_wavevector($cryst, $dims, $demag, q)
    end

    println("\n-- reuse one plan across all q --")
    @btime for q in $qs
        Sunny.precompute_dipole_ewald_at_wavevector($plan, q)
    end
    return nothing
end

run()
