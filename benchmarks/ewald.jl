# Benchmark for the dipole-dipole Ewald builder (src/System/Ewald.jl) on the
# spin-wave / powder workload: a Hamiltonian rebuilt once per wavevector q on a
# single flattened magnetic cell (dims = (1,1,1)).
#
# The physical test case matches PR #505 — a synthetic dipolar magnet on the
# FCC (spacegroup 227) lattice, reshaped to its primitive cell — so the per-q
# timings here are directly comparable to that PR.
#
# Run from the package environment:
#
#     julia --project=. benchmarks/ewald.jl
#
# or interactively: `include("benchmarks/ewald.jl")`.

using Sunny, LinearAlgebra, Random, Printf
using Sunny: ewald_interaction_tensor, EwaldTensorCache, swt_hamiltonian_dipole!, nbands, natoms

# Minimum wall time of `f()` over repeated calls (best-of is the standard robust
# estimator for a deterministic function). Avoids a BenchmarkTools dependency.
function besttime(f; seconds=0.5)
    f()                                     # warm up / compile
    best = Inf
    t0 = time()
    while time() - t0 < seconds
        best = min(best, @elapsed f())
    end
    return best
end

# Quasi-uniform points on a sphere of the given radius (a powder q-shell).
function fibonacci_sphere(n, radius)
    ϕ = (1 + √5) / 2
    return [begin
        z = 1 - 2(i - 0.5) / n
        r = √max(0, 1 - z^2)
        θ = 2π * (i / ϕ)
        radius * Sunny.Vec3(r*cos(θ), r*sin(θ), z)
    end for i in 1:n]
end

# PR #505's synthetic dipolar magnet: FCC (227), s = 7/2, reshaped to the
# primitive cell so the spin-wave calculation runs on a single magnetic cell.
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

# ── Comparison with PR #505 ─────────────────────────────────────────────────
#
# Both approaches precompute the q-independent Ewald pieces once and reuse them
# across all q. This branch stores them in an `EwaldTensorCache` held by
# `sys.ewald` (built once in the `Ewald` constructor); `ewald_interaction_tensor`
# then assembles Aq(q) from the cache, and `swt_hamiltonian_dipole!` calls it
# per q. PR #505 instead caches the local-frame–projected real-space terms as
# fields inside `SpinWaveTheory` (`ewald_q_plan`, `ewald_real_local`), folding
# the local frame (g, √S, R) into the Ewald cache.
#
# The comparison is apples-to-apples: both build the identical FCC-227 s=7/2
# dipolar system reshaped to primitive, sweep the same fibonacci_sphere q-shell,
# and rebuild the same 2L×2L Hamiltonian, timed by the same sequential best-of
# loop. (`swt_hamiltonian_dipole!` has no internal threading, so the process
# thread count is irrelevant; PR #505 additionally threads over q — an orthogonal
# ~Nthreads speedup this per-q kernel benchmark excludes.)
#
# Unlike PR #505, the cache here keeps the Ewald↔SpinWaveTheory boundary intact:
# `ewald_interaction_tensor` returns Aq(q) and SpinWaveTheory owns the local
# frame and does the projection.
function run_ewald_benchmarks(; nq=2000, radius=0.8)
    swt = synthetic_dipolar_swt()
    sys = swt.sys
    na = natoms(sys.crystal)
    L = nbands(swt)
    qs = fibonacci_sphere(nq, radius)

    # (a) The Ewald builder in isolation, reusing the cached q-independent pieces
    # (what src/System/Ewald.jl optimizes).
    cache = sys.ewald.cache
    tewald = besttime() do
        for q in qs
            ewald_interaction_tensor(cache, q)
        end
    end

    # (b) The full spin-wave Hamiltonian rebuild (PR #505's headline metric).
    H = zeros(ComplexF64, 2L, 2L)
    tswt = besttime() do
        for q in qs
            swt_hamiltonian_dipole!(H, swt, q)
        end
    end

    println("PR #505 test case: FCC (227) dipolar, primitive cell, na = $na, L = $L")
    println("Single thread, best-of over $nq wavevectors:\n")
    @printf("  Ewald builder  ewald_interaction_tensor : %6.1f µs/q\n", 1e6*tewald/nq)
    @printf("  Full SWT build swt_hamiltonian_dipole!   : %6.1f µs/q\n", 1e6*tswt/nq)
end

run_ewald_benchmarks()
