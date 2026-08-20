# Tests for the KAExt SU(N) GPU path of SpinWaveTheoryKPM.
#
# Backend: KernelAbstractions.CPU() — runs on CPU for CI portability.
# On CUDA / Metal GPUs the same code paths execute on-device.
#
# Tests:
#   1. SUN matvec regression: batched GPU matvec vs CPU multiply_by_hamiltonian_SUN!
#   2. SUN end-to-end intensities vs Sunny CPU SpinWaveTheoryKPM (SU(2), spin-1/2)
#   3. SUN end-to-end intensities, multi-atom unit cell (SU(3), spin-1)

using Test
using Sunny
using KernelAbstractions
using LinearAlgebra
using Random

backend = KernelAbstractions.CPU()
ka_ext = Base.get_extension(Sunny, :KAExt)
@assert ka_ext !== nothing "KAExt failed to load"

@testset "KA SUN batched matvec vs CPU (SU(2), spin-1/2, single atom)" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=2)], :SUN, seed=1)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure = ssf_perp(sys)
    swt     = SpinWaveTheory(sys; measure)
    Nf      = Sunny.nflavors(swt)
    Na      = Sunny.natoms(swt.sys.crystal)
    L       = Nf * Na
    twoL    = 2L

    # Build device structs
    swt_d    = Sunny.to_device(swt, backend)
    incoming = ka_ext.build_incoming_bond_csr_sun_ka(swt.sys, backend)

    rng = MersenneTwister(42)
    N_chains = 5
    qs = [Sunny.to_reshaped_rlu(swt.sys, Sunny.Vec3(0.1+0.05c, 0.2+0.03c, 0.3+0.07c)) for c in 1:N_chains]
    X  = randn(rng, ComplexF64, N_chains, twoL)

    X_ka = KernelAbstractions.allocate(backend, ComplexF64, N_chains, twoL)
    copyto!(X_ka, X)
    Y_ka = KernelAbstractions.zeros(backend, ComplexF64, N_chains, twoL)
    qs_ka = KernelAbstractions.allocate(backend, Sunny.Vec3, N_chains)
    copyto!(qs_ka, qs)

    ka_ext.mul_dynamical_matrix_batched_ka!(swt_d, incoming, Y_ka, X_ka, qs_ka)
    Y_gpu = Array(Y_ka)

    # CPU reference: use multiply_by_hamiltonian_SUN! per chain
    for c in 1:N_chains
        y_ref = zeros(ComplexF64, 1, twoL)
        Sunny.mul_dynamical_matrix!(swt, y_ref, reshape(X[c, :], 1, :), [qs[c]])
        rel_err = maximum(abs, Y_gpu[c, :] .- vec(y_ref)) / max(maximum(abs, y_ref), 1e-14)
        @test rel_err < 1e-12
    end
end

@testset "KA SUN batched matvec vs CPU (SU(3), spin-1, 2x2x2 supercell)" begin
    # 2×2×2 supercell gives Na=8 atoms after SWT reshaping, exercising multi-atom paths
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :SUN, dims=(2, 2, 2), seed=2)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure = ssf_trace(sys)
    swt     = SpinWaveTheory(sys; measure)
    Nf      = Sunny.nflavors(swt)
    Na      = Sunny.natoms(swt.sys.crystal)
    L       = Nf * Na
    twoL    = 2L

    swt_d    = Sunny.to_device(swt, backend)
    incoming = ka_ext.build_incoming_bond_csr_sun_ka(swt.sys, backend)

    rng = MersenneTwister(7)
    N_chains = 4
    qs = [Sunny.to_reshaped_rlu(swt.sys, Sunny.Vec3(0.15c, 0.25c, 0.1c)) for c in 1:N_chains]
    X  = randn(rng, ComplexF64, N_chains, twoL)

    X_ka = KernelAbstractions.allocate(backend, ComplexF64, N_chains, twoL)
    copyto!(X_ka, X)
    Y_ka = KernelAbstractions.zeros(backend, ComplexF64, N_chains, twoL)
    qs_ka = KernelAbstractions.allocate(backend, Sunny.Vec3, N_chains)
    copyto!(qs_ka, qs)

    ka_ext.mul_dynamical_matrix_batched_ka!(swt_d, incoming, Y_ka, X_ka, qs_ka)
    Y_gpu = Array(Y_ka)

    for c in 1:N_chains
        y_ref = zeros(ComplexF64, 1, twoL)
        Sunny.mul_dynamical_matrix!(swt, y_ref, reshape(X[c, :], 1, :), [qs[c]])
        rel_err = maximum(abs, Y_gpu[c, :] .- vec(y_ref)) / max(maximum(abs, y_ref), 1e-14)
        @test rel_err < 1e-12
    end
end

@testset "KA SUN end-to-end intensities vs CPU (SU(2), spin-1/2)" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=2)], :SUN, seed=3)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure = ssf_perp(sys)
    qs      = [[h, 0.5, 0.5] for h in range(0, 1, length=5)]
    path    = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 5.0, length=8)
    kernel  = gaussian(fwhm=0.3)

    swt_kry  = SpinWaveTheoryKPM(sys; measure, niters=40)
    res_cpu  = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "KA SUN end-to-end intensities vs CPU (SU(3), spin-1, two-atom cell)" begin
    latvecs = lattice_vectors(1.0, 1.0, 2.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0], [0, 0, 0.5]]; types=["A", "B"])
    sys     = System(cryst, [1 => Moment(s=1, g=2), 2 => Moment(s=1, g=2)], :SUN, seed=4)
    set_exchange!(sys, -1.0, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, -0.5, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -0.5, Bond(2, 2, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    minimize_energy!(sys)

    measure  = ssf_trace(sys)
    qs       = [[h, 0.0, 0.0] for h in range(0, 1, length=5)]
    path     = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 8.0, length=10)
    kernel   = gaussian(fwhm=0.3)

    swt_kry  = SpinWaveTheoryKPM(sys; measure, niters=60)
    res_cpu  = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref   = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-6
end

@testset "KA SUN FP32 intensities vs FP64 (SU(2), spin-1/2)" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=2)], :SUN, seed=5)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure   = ssf_perp(sys)
    qs        = [[h, 0.5, 0.5] for h in range(0, 1, length=5)]
    path      = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies  = range(0.1, 5.0, length=8)
    kernel    = gaussian(fwhm=0.3)

    swt_kry   = SpinWaveTheoryKPM(sys; measure, niters=40)
    res_fp64  = Sunny.intensities(Sunny.to_device_batched(swt_kry, backend; precision=Float64),
                                  path; energies, kernel)
    res_fp32  = Sunny.intensities(Sunny.to_device_batched(swt_kry, backend; precision=Float32),
                                  path; energies, kernel)
    maxref    = maximum(abs, res_fp64.data)
    @test maxref > 1e-10
    @test maximum(abs, res_fp32.data .- res_fp64.data) / maxref < 1e-3
end

@testset "KA error: method=:kpm rejected at to_device_batched" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1/2, g=2)], :SUN, seed=6)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    measure  = ssf_perp(sys)
    swt_kry  = SpinWaveTheoryKPM(sys; measure, niters=20, method=:kpm)
    @test_throws ErrorException Sunny.to_device_batched(swt_kry, backend)
end

# ── Tests added during kpm-gpu SU(N) validation ─────────────────────────────
#
# These four testsets cover properties confirmed by the comprehensive benchmark
# validation (benchmarks/sun_comprehensive_validation.jl) that were not covered
# by the testsets above.

@testset "KA SUN odd-niters parity fix (SU(3), s=1, cubic 2×2×2)" begin
    # Fix applied in intensities_lanczos_device_batched_ka!: the FP64 GPU Lanczos
    # path now rounds odd niters up to even (max_iters_global += mod(max_iters_global, 2)),
    # matching Sunny's CPU lanczos() which does the same.  Without the fix, niters=5
    # gave ~27-54% error because CPU silently ran 6 iterations while GPU ran 5.
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :SUN, dims=(2, 2, 2), seed=2)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure  = ssf_trace(sys)
    qs       = [[h, 0.0, 0.0] for h in range(0, 1, length=5)]
    path     = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 8.0, length=8)
    kernel   = gaussian(fwhm=0.3)

    for niters_odd in [5, 7]
        swt_kry = SpinWaveTheoryKPM(sys; measure, niters=niters_odd)
        res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
        maxref  = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        swt_d   = Sunny.to_device_batched(swt_kry, backend)
        res_gpu = Sunny.intensities(swt_d, path; energies, kernel)
        @test maximum(abs, res_gpu.data .- res_cpu.data) / maxref < 1e-5
    end
end

@testset "KA SUN adaptive tol (SU(3), s=1, triangular AFM 3×3)" begin
    latvecs = lattice_vectors(1.0, 1.0, 10.0, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    minimize_energy!(sys)

    measure  = ssf_trace(sys)
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0]], 8)
    energies = range(0.1, 4.0, 8)
    kernel   = lorentzian(fwhm=0.4)

    swt_kry = SpinWaveTheoryKPM(sys; measure, tol=0.05)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref  = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_d   = Sunny.to_device_batched(swt_kry, backend)
    res_gpu = Sunny.intensities(swt_d, path; energies, kernel)
    @test maximum(abs, res_gpu.data .- res_cpu.data) / maxref < 1e-5
end

@testset "KA SUN finite temperature, K→M path (SU(3), s=1, triangular AFM 3×3)" begin
    # kT>0 at q=Γ is numerically unstable: the Goldstone mode has a Ritz eigenvalue
    # ~±1e-8, and thermal_prefactor sensitivity kT/λ² ≈ 5×10¹² at (λ=1e-8, kT=0.5)
    # amplifies a 2×10⁻¹³ CPU/GPU FP difference to ~1000 thermal-weight error → 16%.
    # A K→M path (no Γ) avoids this; CPU and GPU agree to machine precision there.
    latvecs = lattice_vectors(1.0, 1.0, 10.0, 90, 90, 120)
    cryst   = Crystal(latvecs, [[0, 0, 0]])
    sys     = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims=(3, 3, 1), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    minimize_energy!(sys)

    measure  = ssf_trace(sys)
    path     = q_space_path(cryst, [[1/3,1/3,0],[1/2,0,0],[1/3,1/3,0]], 8)
    energies = range(0.1, 4.0, 8)
    kernel   = lorentzian(fwhm=0.4)

    for kT in [0.5, 2.0]
        swt_kry = SpinWaveTheoryKPM(sys; measure, niters=10)
        res_cpu = Sunny.intensities(swt_kry, path; energies, kernel, kT)
        maxref  = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        swt_d   = Sunny.to_device_batched(swt_kry, backend)
        res_gpu = Sunny.intensities(swt_d, path; energies, kernel, kT)
        @test maximum(abs, res_gpu.data .- res_cpu.data) / maxref < 1e-5
    end
end

@testset "KA SUN scale independence niters=10 (SU(3), s=1, triangular AFM)" begin
    # GPU kernel accuracy at niters=10 (~1e-5 for SU(N) bandwidth ~24) must be
    # independent of system size Na.  Tests Na=9 (3×3×1) and Na=36 (6×6×1).
    latvecs  = lattice_vectors(1.0, 1.0, 10.0, 90, 90, 120)
    cryst    = Crystal(latvecs, [[0, 0, 0]])
    path     = q_space_path(cryst, [[0,0,0],[1/3,1/3,0],[1/2,0,0]], 8)
    energies = range(0.1, 4.0, 8)
    kernel   = lorentzian(fwhm=0.4)

    for dims in [(3, 3, 1), (6, 6, 1)]
        sys = System(cryst, [1 => Moment(s=1, g=2)], :SUN; dims, seed=0)
        set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
        polarize_spins!(sys, [0, 0, 1])
        minimize_energy!(sys)

        measure = ssf_trace(sys)
        swt_kry = SpinWaveTheoryKPM(sys; measure, niters=10)
        res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
        maxref  = maximum(abs, res_cpu.data)
        @test maxref > 1e-10

        swt_d   = Sunny.to_device_batched(swt_kry, backend)
        res_gpu = Sunny.intensities(swt_d, path; energies, kernel)
        @test maximum(abs, res_gpu.data .- res_cpu.data) / maxref < 1e-5
    end
end
