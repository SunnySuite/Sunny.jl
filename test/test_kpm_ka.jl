# Tests for the KAExt FP64 GPU path of SpinWaveTheoryKPM.
#
# Backend:  KernelAbstractions.CPU() — runs the FP64 driver on CPU for CI
#           portability.  On CUDA / AMD GPUs the same code paths execute
#           on-device at native Float64 precision.
#
# All GPU paths use Lanczos (method=:kpm is rejected at to_device_batched time).
# (FP32 dispatch is exercised in test/test_kpm_ka_fp32.jl.)
#
# Reference for every accuracy test is Sunny's unmodified
# `intensities(::SpinWaveTheoryKPM)` on the CPU.

using Test
using Sunny
using KernelAbstractions
using LinearAlgebra
using Random

backend = KernelAbstractions.CPU()
ka_ext = Base.get_extension(Sunny, :KAExt)
@assert ka_ext !== nothing "KAExt failed to load"

@testset "KA matvec per-q vs CPU" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure = ssf_perp(sys)
    swt = SpinWaveTheory(sys; measure)
    L = Sunny.natoms(swt.sys.crystal)
    twoL = 2L
    swt_d = Sunny.to_device(swt, backend)

    rng = MersenneTwister(42)
    q = Sunny.Vec3(0.2, 0.3, 0.1)
    q_reshaped = Sunny.to_reshaped_rlu(swt.sys, q)
    x = randn(rng, ComplexF64, twoL)

    y_cpu = zeros(ComplexF64, twoL)
    Sunny.mul_dynamical_matrix!(swt, reshape(y_cpu, 1, :), reshape(x, 1, :), [q_reshaped])

    x_ka = KernelAbstractions.allocate(backend, ComplexF64, 1, twoL)
    copyto!(x_ka, reshape(x, 1, :))
    y_ka = KernelAbstractions.zeros(backend, ComplexF64, 1, twoL)
    qs_ka = KernelAbstractions.allocate(backend, Sunny.Vec3, 1)
    copyto!(qs_ka, [q_reshaped])
    Sunny.mul_dynamical_matrix!(swt_d, y_ka, x_ka, qs_ka)

    @test maximum(abs, Array(y_ka)[1, :] .- y_cpu) / maximum(abs, y_cpu) < 1e-13
end

@testset "KA batched matvec vs per-q" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])

    measure = ssf_perp(sys)
    swt = SpinWaveTheory(sys; measure)
    L = Sunny.natoms(swt.sys.crystal)
    twoL = 2L
    swt_d = Sunny.to_device(swt, backend)
    incoming = ka_ext.build_incoming_bond_csr_ka(swt.sys, backend)

    rng = MersenneTwister(7)
    N_chains = 4
    qs = [Sunny.to_reshaped_rlu(swt.sys, Sunny.Vec3(0.1*c, 0.2*c, 0.3*c)) for c in 1:N_chains]
    X = randn(rng, ComplexF64, N_chains, twoL)

    X_ka = KernelAbstractions.allocate(backend, ComplexF64, N_chains, twoL); copyto!(X_ka, X)
    Y_ka = KernelAbstractions.zeros(backend, ComplexF64, N_chains, twoL)
    qs_ka = KernelAbstractions.allocate(backend, Sunny.Vec3, N_chains); copyto!(qs_ka, qs)
    ka_ext.mul_dynamical_matrix_batched_ka!(swt_d, incoming, Y_ka, X_ka, qs_ka)
    Y_batched = Array(Y_ka)

    for c in 1:N_chains
        y_ref = zeros(ComplexF64, 1, twoL)
        Sunny.mul_dynamical_matrix!(swt, y_ref, reshape(X[c, :], 1, :), [qs[c]])
        @test maximum(abs, Y_batched[c, :] .- vec(y_ref)) / maximum(abs, y_ref) < 1e-13
    end
end

@testset "KA end-to-end intensities vs Sunny CPU (cubic FM primitive)" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)

    qs = [[h, 0.5, 0.5] for h in range(0, 1, length=5)]
    path = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 8.0, length=8)
    kernel = gaussian(fwhm=0.3)

    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=40)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_perq  = Sunny.to_device(swt_kry, backend)
    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_perq  = Sunny.intensities(swt_perq,  path; energies, kernel)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_perq.data  .- res_cpu.data) / maxref < 1e-10
    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "KA end-to-end intensities vs Sunny CPU (cubic FM 2x2x2)" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)

    qs = [[h, 0.5, 0.5] for h in range(0, 1, length=5)]
    path = q_space_path(Sunny.orig_crystal(sys), qs, 5)
    energies = range(0.1, 8.0, length=10)
    kernel = gaussian(fwhm=0.3)

    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=60)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_perq  = Sunny.to_device(swt_kry, backend)
    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_perq  = Sunny.intensities(swt_perq,  path; energies, kernel)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_perq.data  .- res_cpu.data) / maxref < 1e-10
    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "KA with biquadratic coupling" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]); biquad=-0.1)
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]); biquad=-0.1)
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]); biquad=-0.1)
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)

    qs = [[h, 0.5, 0.5] for h in range(0, 1, length=4)]
    path = q_space_path(Sunny.orig_crystal(sys), qs, 4)
    energies = range(0.1, 8.0, length=8)
    kernel = gaussian(fwhm=0.3)

    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=60)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_perq  = Sunny.to_device(swt_kry, backend)
    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_perq  = Sunny.intensities(swt_perq,  path; energies, kernel)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_perq.data  .- res_cpu.data) / maxref < 1e-10
    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "KA with different q-path direction" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)

    qs = [[h, h, h] for h in range(0.1, 0.9, length=4)]
    path = q_space_path(Sunny.orig_crystal(sys), qs, 4)
    energies = range(0.1, 8.0, length=8)
    kernel = gaussian(fwhm=0.3)

    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=60)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_perq  = Sunny.to_device(swt_kry, backend)
    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_perq  = Sunny.intensities(swt_perq,  path; energies, kernel)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_perq.data  .- res_cpu.data) / maxref < 1e-10
    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "KA with ssf_trace measure" begin
    latvecs = lattice_vectors(1, 1, 1, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 1, 0]))
    set_exchange!(sys, -1.0, Bond(1, 1, [0, 0, 1]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_trace(sys)

    qs = [[h, 0.5, 0.5] for h in range(0, 1, length=4)]
    path = q_space_path(Sunny.orig_crystal(sys), qs, 4)
    energies = range(0.1, 8.0, length=8)
    kernel = gaussian(fwhm=0.3)

    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=60)
    res_cpu = Sunny.intensities(swt_kry, path; energies, kernel)
    maxref = maximum(abs, res_cpu.data)
    @test maxref > 1e-10

    swt_perq  = Sunny.to_device(swt_kry, backend)
    swt_batch = Sunny.to_device_batched(swt_kry, backend)
    res_perq  = Sunny.intensities(swt_perq,  path; energies, kernel)
    res_batch = Sunny.intensities(swt_batch, path; energies, kernel)

    @test maximum(abs, res_perq.data  .- res_cpu.data) / maxref < 1e-10
    @test maximum(abs, res_batch.data .- res_cpu.data) / maxref < 1e-10
end

@testset "method=:kpm rejected at to_device_batched" begin
    latvecs = lattice_vectors(1.0, 1.0, 1.0, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, dims=(2,2,2), seed=0)
    set_exchange!(sys, -1.0, Bond(1, 1, [1, 0, 0]))
    polarize_spins!(sys, [0, 0, 1])
    measure = ssf_perp(sys)
    swt_kry = SpinWaveTheoryKPM(sys; measure, niters=60, method=:kpm)
    @test_throws ErrorException Sunny.to_device_batched(swt_kry, backend)
end
