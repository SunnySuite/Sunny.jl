# Tests for the KAExt FP32 GPU path of SpinWaveTheoryKPM.
#
# Backend:  KernelAbstractions.CPU() — runs the FP32 driver on CPU for CI
#           portability.  On Metal, the same code paths execute on Apple
#           Silicon GPUs at native Float32 precision; see the benchmark
#           file `benchmarks/mac_metal_lanczos_benchmark.jl` for hardware
#           speedup numbers.
#
# All GPU paths use Lanczos; method=:kpm is rejected at to_device_batched time.
#
# Reference for every accuracy test is Sunny's unmodified
# `intensities(::SpinWaveTheoryKPM)` on the CPU.
#
# Regularization note:  FP32 GPU Lanczos uses `regularization=1e-4` for
# tests that exercise systems with spectral degeneracies (homogeneous
# higher-spin or biquadratic-only couplings).  The Sunny default of 1e-8
# is FP64-safe but, in FP32, those degeneracies amplify rounding error
# into ghost eigenvalues.  Bumping the regularization opens a small gap
# that puts the spectrum above the FP32 noise floor.  See:
# `SpinWaveTheoryKPM` docstring and the benchmark file for guidance.

using Test
using Sunny
using KernelAbstractions
using LinearAlgebra
using Random

backend = KernelAbstractions.CPU()
ka_ext = Base.get_extension(Sunny, :KAExt)
@assert ka_ext !== nothing "KAExt failed to load"

const FP32_RTOL  = 1e-3   # FP32 production tolerance (per validation matrix)
const FP32_REG   = 1e-4   # Recommended regularization for FP32 (see notes above)
const SEED       = 42

# Build a SpinWaveTheoryKPM with the FP32-recommended regularization.
mkswt(sys, measure; method=:lanczos, kwargs...) =
    SpinWaveTheoryKPM(sys; measure, regularization=FP32_REG, method, kwargs...)


# ---------------------------------------------------------------------------
# API: dispatch and rejection of unsupported precisions
# ---------------------------------------------------------------------------

@testset "KA: to_device_batched precision dispatch" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :dipole; seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    swt = SpinWaveTheoryKPM(sys; measure=ssf_perp(sys), tol=0.05, regularization=FP32_REG)

    swt_fp64_default  = Sunny.to_device_batched(swt, backend)
    swt_fp64_explicit = Sunny.to_device_batched(swt, backend; precision=Float64)
    swt_fp32          = Sunny.to_device_batched(swt, backend; precision=Float32)

    @test swt_fp64_default.precision  === Float64
    @test swt_fp64_explicit.precision === Float64
    @test swt_fp32.precision          === Float32

    @test_throws ErrorException Sunny.to_device_batched(swt, backend; precision=Float16)
    @test_throws ErrorException Sunny.to_device_batched(swt, backend; precision=BigFloat)
end

@testset "KA: rejects unsupported method" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :dipole; seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    @test_throws ErrorException Sunny.to_device_batched(
        SpinWaveTheoryKPM(sys; measure=ssf_perp(sys), tol=0.05, method=:bogus),
        backend; precision=Float32,
    )
end


# ---------------------------------------------------------------------------
# FP32 Lanczos accuracy vs CPU FP64 Lanczos (apples-to-apples; same algorithm)
# ---------------------------------------------------------------------------

@testset "FP32 Lanczos: cubic FM s=1, fixed niters" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :dipole; dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)

    qs       = [[h, 0.5, 0.5] for h in range(0, 1, length=5)]
    path     = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 8.0, length=10)
    kernel   = gaussian(fwhm=0.3)

    swt = mkswt(sys, measure; niters=40)  # even, well below natural-convergence cap
    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: cubic s=2 (regularization=1e-4 stabilizes)" begin
    # Higher-spin homogeneous system — fails FP32 with default reg=1e-8 due
    # to spectral degeneracies, passes with reg=1e-4 (mentor's recommendation).
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=2, g=2)], :dipole; dims=(4, 4, 1), seed=0)
    set_exchange!(sys, 1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(sys)

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/2,0,0],[1/2,1/2,0],[0,0,0]], 30)
    energies = range(0, 20, 50)
    kernel   = lorentzian(fwhm=1.0)

    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: biquadratic s=1 (regularization=1e-4 stabilizes)" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=1)], :dipole; dims=(2, 2, 1), seed=0)
    set_exchange!(sys, 0.7, Bond(1, 1, [1, 0, 0]); biquad=0.3)
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/2,0,0],[1/2,1/2,0],[0,0,0]], 30)
    energies = range(0, 8, 50)
    kernel   = lorentzian(fwhm=0.5)

    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: disordered triangular AFM (Ex09 30x30) tol-based" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (10, 10, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0],[0,0,0]], 60)
    energies = range(0.0, 3.0, 60)
    kernel   = lorentzian(fwhm=0.4)

    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: external field opens spectral gap" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    set_field!(sys_in, [0, 0, 7.5])
    randomize_spins!(sys_in); minimize_energy!(sys_in, maxiters=2_000)

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0],[0,0,0]], 30)
    energies = range(0.0, 9.0, 30)
    kernel   = lorentzian(fwhm=0.4)

    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end


# ---------------------------------------------------------------------------
# FP32 Lanczos capabilities not available in FP32 Chebyshev:
#   - finite temperature kT > 0
#   - odd user-niters (rounded up to even, matching CPU Lanczos)
# ---------------------------------------------------------------------------

@testset "FP32 Lanczos: finite kT (kT=0.5)" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0]], 20)
    energies = range(0.0, 3.0, 30)
    kernel   = lorentzian(fwhm=0.4)

    Random.seed!(SEED)
    res_cpu = Sunny.intensities(swt, path; energies, kernel, kT=0.5)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel, kT=0.5)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: odd user-niters round to even (parity fix)" begin
    # Without the parity fix, odd user-specified niters produced
    # catastrophic FP32 errors (rel err ~0.9).  The driver now rounds
    # niters up to the nearest even integer, matching CPU Lanczos's
    # internal `niters += mod(niters, 2)` adjustment.
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    measure  = ssf_perp(sys_in)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0]], 20)
    energies = range(0.0, 3.0, 30)
    kernel   = lorentzian(fwhm=0.4)

    for n_user in (25, 35)   # both odd; would have failed pre-parity-fix
        swt = mkswt(sys_in, measure; niters=n_user)
        Random.seed!(SEED)
        res_cpu = Sunny.intensities(swt, path; energies, kernel)
        maxref  = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
        Random.seed!(SEED)
        res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

        @test all(isfinite, res_fp32.data)
        @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
    end
end


# ---------------------------------------------------------------------------
# Composition: the FP32 GPU wrapper drops into Sunny.powder_average with
# zero source changes.  Both paths go through the identical
# `powder_average` driver (same seed → identical Fibonacci sampling and
# identical random shell rotations); only the precision of the inner
# `intensities(swt, qs; ...)` call differs.  Confirms that research
# workflows built on top of `powder_average` transparently inherit the
# GPU speedup.
# ---------------------------------------------------------------------------

@testset "FP32 Lanczos: powder_average composition" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (3, 3, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    swt      = mkswt(sys_in, ssf_perp(sys_in); tol=0.05)
    energies = range(0.0, 3.0, 10)
    kernel   = lorentzian(fwhm=0.4)
    radii    = range(0.5, 2.5, 3)     # 3 radial shells
    n_fib    = 30                      # 30 Fibonacci points per shell (90 q total)

    Random.seed!(SEED)
    res_cpu = powder_average(cryst, radii, n_fib) do qs
        intensities(swt, qs; energies, kernel)
    end
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = powder_average(cryst, radii, n_fib) do qs
        Sunny.intensities(swt_fp32, qs; energies, kernel)
    end

    @test size(res_fp32.data) == size(res_cpu.data)
    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end



# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

@testset "FP32 Lanczos: 3-point path" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    measure  = ssf_perp(sys_in)
    swt      = mkswt(sys_in, measure; tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0]], 3)  # short path
    energies = range(0.0, 3.0, 30)
    kernel   = lorentzian(fwhm=0.4)

    Random.seed!(SEED)
    res_cpu  = Sunny.intensities(swt, path; energies, kernel)
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos: ssf_trace observable" begin
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=1)], :dipole; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, +1.0, Bond(1, 1, [1, 0, 0]))
    randomize_spins!(sys); minimize_energy!(sys)
    sys_in  = to_inhomogeneous(repeat_periodically(sys, (5, 5, 1)))
    for (s1, s2, off) in symmetry_equivalent_bonds(sys_in, Bond(1, 1, [1, 0, 0]))
        set_exchange_at!(sys_in, 1.0 + randn(sys_in.rng)/3, s1, s2; offset=off)
    end
    minimize_energy!(sys_in, maxiters=2_000)

    swt      = mkswt(sys_in, ssf_trace(sys_in); tol=0.05)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0]], 20)
    energies = range(0.0, 3.0, 30)
    kernel   = lorentzian(fwhm=0.4)

    Random.seed!(SEED)
    res_cpu  = Sunny.intensities(swt, path; energies, kernel)
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, path; energies, kernel)

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end


# ---------------------------------------------------------------------------
# SU(N) FP32 Lanczos — Metal production path
#
# Metal cannot execute Float64 arithmetic inside kernels, so FP32
# (precision=Float32) is the only viable GPU path on Apple Silicon.
# All testsets below validate the FP32 SUN kernel against the CPU FP64
# reference (same swt object, same niters/tol) at FP32_RTOL = 1e-3.
#
# System: triangular FM, s=1, SU(3) (Nf=2), dims=(3,3,1) → Na=9, twoL=36.
# The 3×3 supercell is the minimal cell that accommodates the √3×√3
# triangular ground state; minimize_energy! is required to reach it.
# ---------------------------------------------------------------------------

# Shared SU(N) fixture used by all SUN FP32 testsets below.
const SUN_CRYST   = Crystal(lattice_vectors(1, 1, 10, 90, 90, 120), [[0, 0, 0]])
const SUN_SYS     = let
    s = System(SUN_CRYST, [1 => Moment(s=1, g=2)], :SUN; dims=(3, 3, 1), seed=0)
    set_exchange!(s, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(s, [0, 0, 1])
    minimize_energy!(s)
    s
end
const SUN_PATH    = q_space_path(SUN_CRYST, [[0,0,0],[1/3,1/3,0],[1/2,0,0],[0,0,0]], 15)
const SUN_PATH_KT = q_space_path(SUN_CRYST, [[1/3,1/3,0],[1/2,0,0],[1/3,1/3,0]], 10)
const SUN_ENERGIES = range(0.0, 4.0, 10)

@testset "FP32 Lanczos SUN: basic accuracy (SU(3), s=1, 3×3×1, niters=10)" begin
    # Core correctness: FP32 SUN kernel vs CPU FP64 reference.
    # niters=10 is well below the Paige ghost region for twoL=36.
    measure  = ssf_trace(SUN_SYS)
    swt      = mkswt(SUN_SYS, measure; niters=10)

    Random.seed!(SEED)
    res_cpu  = Sunny.intensities(swt, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos SUN: ssf_perp measure" begin
    measure  = ssf_perp(SUN_SYS)
    swt      = mkswt(SUN_SYS, measure; niters=10)

    Random.seed!(SEED)
    res_cpu  = Sunny.intensities(swt, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos SUN: broadening kernels (gaussian, lorentzian wide/narrow)" begin
    measure  = ssf_trace(SUN_SYS)
    swt      = mkswt(SUN_SYS, measure; niters=10)
    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)

    for kern in [gaussian(fwhm=0.4), lorentzian(fwhm=1.5), lorentzian(fwhm=0.1)]
        Random.seed!(SEED)
        res_cpu  = Sunny.intensities(swt, SUN_PATH; energies=SUN_ENERGIES, kernel=kern)
        maxref   = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        Random.seed!(SEED)
        res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH; energies=SUN_ENERGIES, kernel=kern)

        @test all(isfinite, res_fp32.data)
        @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
    end
end

@testset "FP32 Lanczos SUN: adaptive tol=0.05" begin
    measure  = ssf_trace(SUN_SYS)
    swt      = mkswt(SUN_SYS, measure; tol=0.05)

    Random.seed!(SEED)
    res_cpu  = Sunny.intensities(swt, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
    Random.seed!(SEED)
    res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))

    @test all(isfinite, res_fp32.data)
    @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
end

@testset "FP32 Lanczos SUN: finite kT, K→M path (no Goldstone singularity)" begin
    # q=Γ has a Goldstone mode with Ritz eigenvalue ~±1e-8; thermal_prefactor
    # sensitivity kT/λ² ≈ 5×10¹² amplifies even FP32 noise there.
    # K→M→K avoids Γ; all off-Γ q-points agree within FP32_RTOL.
    measure  = ssf_trace(SUN_SYS)
    swt      = mkswt(SUN_SYS, measure; niters=10)
    swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)

    for kT in [0.5, 2.0]
        Random.seed!(SEED)
        res_cpu  = Sunny.intensities(swt, SUN_PATH_KT; energies=SUN_ENERGIES,
                                     kernel=lorentzian(fwhm=0.4), kT)
        maxref   = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        Random.seed!(SEED)
        res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH_KT; energies=SUN_ENERGIES,
                                     kernel=lorentzian(fwhm=0.4), kT)

        @test all(isfinite, res_fp32.data)
        @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
    end
end

@testset "FP32 Lanczos SUN: odd user-niters parity fix" begin
    # The FP32 driver rounds odd niters up to even (same fix as FP64 path).
    # Without the fix, niters=5 produced ~27-54% error for SUN mode.
    measure = ssf_trace(SUN_SYS)

    for niters_odd in [5, 7]
        swt      = mkswt(SUN_SYS, measure; niters=niters_odd)

        Random.seed!(SEED)
        res_cpu  = Sunny.intensities(swt, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))
        maxref   = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        swt_fp32 = Sunny.to_device_batched(swt, backend; precision=Float32)
        Random.seed!(SEED)
        res_fp32 = Sunny.intensities(swt_fp32, SUN_PATH; energies=SUN_ENERGIES, kernel=lorentzian(fwhm=0.4))

        @test all(isfinite, res_fp32.data)
        @test maximum(abs, res_fp32.data .- res_cpu.data) / maxref < FP32_RTOL
    end
end
