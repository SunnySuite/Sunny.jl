import LinearAlgebra: BlasInt

# Storage for `bogoliubov!`, holding the transformation `T` and energies `ε`
# that it produces, together with the LAPACK scratch space that `eigen!` would
# otherwise allocate on every call.
struct BogoliubovWorkspace
    T     :: Matrix{ComplexF64}   # Transformation, as documented by `excitations!`
    ε     :: Vector{Float64}      # Signed energies, in the same grouping
    Tscr  :: Matrix{ComplexF64}   # Scratch for the reordering below
    εscr  :: Vector{Float64}
    work  :: Vector{ComplexF64}   # LAPACK scratch for zhegvd
    rwork :: Vector{Float64}
    iwork :: Vector{BlasInt}
    info  :: Base.RefValue{BlasInt}
end

function BogoliubovWorkspace(L::Int)
    n = 2L
    # Minimum scratch sizes documented for ZHEGVD with jobz = 'V' and n > 1
    # (always true here, since n = 2L).
    return BogoliubovWorkspace(zeros(ComplexF64, n, n), zeros(n), zeros(ComplexF64, n, n), zeros(n),
                               zeros(ComplexF64, 2n + n^2), zeros(1 + 5n + 2n^2), zeros(BlasInt, 3 + 5n),
                               Ref(BlasInt(0)))
end

# Bogoliubov transformation that diagonalizes a quadratic bosonic Hamiltonian,
# allowing for anomalous terms. The general procedure derives from Colpa,
# Physica A, 93A, 327-353 (1978). Overwrites data in H. The returned energies
# alias `ws.ε`, and so are valid only until the next call on the same workspace.
function bogoliubov!(ws::BogoliubovWorkspace, H::Matrix{ComplexF64})
    (; T, ε, Tscr, εscr, work, rwork, iwork, info) = ws
    n = size(T, 1)
    L = div(n, 2)
    @assert size(H) == (n, n)
    # `eigen!` rejects a non-finite matrix before reaching LAPACK, which is not
    # guaranteed to terminate on one. Retain that check.
    all(isfinite, H) || throw(ArgumentError("Hamiltonian contains non-finite entries"))

    # Initialize T to the para-unitary identity Ĩ = diagm([ones(L), -ones(L)])
    T .= 0
    for i in 1:L
        T[i, i] = 1
        T[i+L, i+L] = -1
    end

    # Solve generalized eigenvalue problem, Ĩ t = λ H t, for columns t of T. This is
    # `eigen!(Hermitian(T), Hermitian(H))` with every allocation hoisted into `ws`.
    ccall((BLAS.@blasfunc(zhegvd_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
           Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
           Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
           Clong, Clong),
          1, 'V', 'U', n, T, n, H, n, ε, work, length(work), rwork,
          length(rwork), iwork, length(iwork), info, 1, 1)
    info[] < 0 && throw(ArgumentError("Invalid argument #$(-info[]) to LAPACK zhegvd"))
    info[] > 0 && throw(PosDefException(info[]))

    # By Sylvester's theorem, "inertia" (sign signature) is invariant under a
    # congruence transform Ĩ → √H Ĩ √H, so exactly L of the λ are negative. LAPACK
    # returns them in ascending order, whereas they are wanted with the positive values
    # first and otherwise ascending in absolute value. That reordering is therefore the
    # fixed permutation [L+1:2L; L:-1:1], and no sort is required. It is not a product
    # of disjoint transpositions, so both T and ε are permuted out of scratch copies.
    #
    # Degenerate λ are the one case where this differs from sorting: `eigen!` breaks ties
    # with an unstable QuickSort, so the order within a degenerate block used to be
    # arbitrary, whereas it is now inherited from LAPACK. Either choice is a valid
    # eigenbasis, differing by a rotation within the block, and the permutation below is
    # at least reproducible.
    @assert ε[L] < 0 < ε[L+1]
    copyto!(Tscr, T)
    copyto!(εscr, ε)
    for j in 1:L
        for i in 1:n
            T[i, j]   = Tscr[i, L+j]
            T[i, L+j] = Tscr[i, L+1-j]
        end
        ε[j]   = εscr[L+j]
        ε[L+j] = εscr[L+1-j]
    end

    # Normalize columns of T so that para-unitarity holds, T† Ĩ T = Ĩ.
    for j in 1:n
        c = 1 / sqrt(abs(ε[j]))
        for i in 1:n
            T[i, j] *= c
        end
    end

    # Inverse of λ are eigenvalues of Ĩ H, or equivalently, of √H Ĩ √H. The first L
    # elements are positive, while the next L elements are negative. Their absolute
    # values are excitation energies for the wavevectors q and -q, respectively.
    @. ε = 1 / ε

    return ε
end

# Variant for callers outside a hot loop, which allocates a workspace per call. Unlike
# the method above, the returned energies are freshly allocated, and so remain valid.
function bogoliubov!(T::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    L = div(size(H, 1), 2)
    @assert size(T) == size(H) == (2L, 2L)
    ws = BogoliubovWorkspace(L)
    bogoliubov!(ws, H)
    copyto!(T, ws.T)
    return ws.ε
end


# Returns |1 + nB(ω)| where nB(ω) = 1 / (exp(βω) - 1) is the Bose function.
# Equivalent to |1 / expm1(-βω)| where expm1(x) = e^x-1.
function thermal_prefactor(ω; kT)
    @assert kT >= 0
    iszero(ω) && return Inf
    return abs(1 / expm1(-ω/kT))
end


"""
    excitations!(T, tmp, swt::SpinWaveTheory, q)

Given a wavevector `q`, solves for the matrix `T` representing quasi-particle
excitations, and returns a list of quasi-particle energies. Both `T` and `tmp`
must be supplied as ``2L×2L`` complex matrices, where ``L`` is the number of
bands for a single ``𝐪`` value.

The columns of `T` are understood to be contracted with the Holstein-Primakoff
bosons ``[𝐛_𝐪, 𝐛_{-𝐪}^†]``. The first ``L`` columns provide the eigenvectors
of the quadratic Hamiltonian for the wavevector ``𝐪``. The next ``L`` columns
of `T` describe eigenvectors for ``-𝐪``. The return value is a vector with
similar grouping: the first ``L`` values are energies for ``𝐪``, and the next
``L`` values are the _negation_ of energies for ``-𝐪``.

    excitations!(T, tmp, swt::SpinWaveTheorySpiral, q; branch)

Calculations on a [`SpinWaveTheorySpiral`](@ref) additionally require a `branch`
index. The possible branches ``(1, 2, 3)`` correspond to scattering processes
``𝐪 - 𝐤, 𝐪, 𝐪 + 𝐤`` respectively, where ``𝐤`` is the ordering wavevector.
Each branch will contribute ``L`` excitations, where ``L`` is the number of
spins in the magnetic cell. This yields a total of ``3L`` excitations for a
given momentum transfer ``𝐪``.
"""
function excitations!(T, tmp, swt::SpinWaveTheory, q)
    L = nbands(swt)
    size(T) == size(tmp) == (2L, 2L) || error("Arguments T and tmp must be $(2L)×$(2L) matrices")
    ws = BogoliubovWorkspace(L)
    energies = excitations!(ws, tmp, swt, q)
    copyto!(T, ws.T)
    return energies
end

# Allocation-free variant that stores the transformation in `ws.T`. The returned
# energies alias `ws.ε`.
function excitations!(ws::BogoliubovWorkspace, H, swt::SpinWaveTheory, q)
    q_reshaped = to_reshaped_rlu(swt.sys, q)
    dynamical_matrix!(H, swt, q_reshaped)

    try
        return bogoliubov!(ws, H)
    catch err
        if err isa PosDefException
            rethrow(InstabilityError("Not an energy-minimum; wavevector q = $(vec3_to_string(q)) unstable."))
        else
            rethrow(err)
        end
    end
end

"""
    excitations(swt::SpinWaveTheory, q)
    excitations(swt::SpinWaveTheorySpiral, q; branch)

Returns a pair `(energies, T)` providing the excitation energies and
eigenvectors. Prefer [`excitations!`](@ref) for performance, which avoids matrix
allocations. See the documentation of [`excitations!`](@ref) for more details.
"""
function excitations(swt::SpinWaveTheory, q)
    L = nbands(swt)
    T = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    energies = excitations!(T, copy(H), swt, q)
    return (energies, T)
end

"""
    dispersion(swt::SpinWaveTheory, qpts)

Given a list of wavevectors `qpts` in reciprocal lattice units (RLU), returns
excitation energies for each band. The return value `ret` is 2D array, and
should be indexed as `ret[band_index, q_index]`.
"""
function dispersion(swt::SpinWaveTheory, qpts)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    disp = zeros(L, length(qpts.qs))
    for (iq, q) in enumerate(qpts.qs)
        view(disp, :, iq) .= view(excitations(swt, q)[1], 1:L)
    end
    return reshape(disp, L, size(qpts.qs)...)
end

"""
    intensities_bands(swt::SpinWaveTheory, qpts; kT=0, threaded=false)

Calculate spin wave excitation bands for a set of ``𝐪``-points in reciprocal
space. This calculation is analogous to [`intensities`](@ref), but does not
perform line broadening of the bands. Use [`load_blas_for_threading`](@ref) and
set `threaded=true` to parallelize the calculation using Julia threads.
"""
function intensities_bands(swt::SpinWaveTheory, qpts; kT=0, with_negative=false, threaded=false)
    (; sys, measure) = swt
    num_observables(measure) == 0 && error("No observables! Construct SpinWaveTheory with a `measure` argument.")
    with_negative && error("Option `with_negative=true` not yet supported.")

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    # Number of (magnetic) atoms in magnetic cell
    @assert sys.dims == (1,1,1)
    Na = nsites(sys)
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Temporary storage for pair correlations
    Nobs = num_observables(measure)
    Ncorr = num_correlations(measure)

    disp = zeros(Float64, L, Nq)
    data = zeros(eltype(measure), L, Nq)

    newbuf() = (
        ws   = BogoliubovWorkspace(L),
        H    = zeros(ComplexF64, 2L, 2L),
        u    = zeros(ComplexF64, 2L, Nobs),
        Avec = zeros(ComplexF64, Nobs),
        corr = zeros(ComplexF64, Ncorr),
    )

    function calc_iq!(buf, iq)
        (; ws, H, u, Avec, corr) = buf
        T = ws.T
        q = qpts.qs[iq]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q
        view(disp, :, iq) .= view(excitations!(ws, H, swt, q), 1:L)

        # Linearized observables Â_μ(q) in Holstein-Primakoff bosons
        set_swt_observable_vectors!(u, swt, q_reshaped, q_global)

        # Fill `intensity` array
        for band in 1:L
            for μ in 1:Nobs
                # Left matrix amplitude ⟨0|Â_μ(q)†|n⟩
                Avec[μ] = dot(view(u, :, μ), view(T, :, band))
            end
            map!(corr, measure.corr_pairs) do (μ, ν)
                # Pair correlations ⟨0|Â_μ(q)†|n⟩⟨n|Â_ν(q)|0⟩
                Avec[μ] * conj(Avec[ν]) / Ncells
            end
            data[band, iq] = thermal_prefactor(disp[band, iq]; kT) * measure.combiner(q_global, corr)
        end
    end

    foreach_chunked(calc_iq!, newbuf, 1:Nq; threaded, warn_blas=true)

    # Reassigning `disp` and `data` would box these captured variables
    disp_reshaped = reshape(disp, L, size(qpts.qs)...)
    data_reshaped = reshape(data, L, size(qpts.qs)...)

    return BandIntensities(cryst, qpts, disp_reshaped, data_reshaped)
end

"""
    intensities!(data, swt::SpinWaveTheory, qpts; energies, kernel, kT=0, threaded=false)
    intensities!(data, sc::SampledCorrelations, qpts; energies, kernel=nothing, kT=0)

Like [`intensities`](@ref), but makes use of storage space `data` to avoid
allocation costs.
"""
function intensities!(data, swt::AbstractSpinWaveTheory, qpts; energies, kernel::AbstractBroadening, kT=0, threaded=false)
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), size(qpts.qs)...)
    bands = intensities_bands(swt, qpts; kT, threaded)
    @assert eltype(bands) == eltype(data)
    broaden!(data, bands; energies, kernel, threaded)
    return Intensities(bands.crystal, bands.qpts, collect(Float64, energies), data)
end

"""
    intensities(swt::SpinWaveTheory, qpts; energies, kernel, kT=0, threaded=false)
    intensities(sc::SampledCorrelations, qpts; energies, kernel=nothing, kT)

Calculates dynamical pair correlation intensities for a set of ``𝐪``-points in
reciprocal space.

Linear spin wave theory calculations are performed with an instance of
[`SpinWaveTheory`](@ref). The alternative [`SpinWaveTheorySpiral`](@ref) allows
to study generalized spiral orders with a single, incommensurate-``𝐤`` ordering
wavevector. Another alternative [`SpinWaveTheoryKPM`](@ref) is favorable for
calculations on large magnetic cells, and allows to study systems with disorder.
An optional nonzero temperature `kT` will scale intensities by the quantum
thermal occupation factor ``|1 + n_B(ω)|`` where ``n_B(ω) = 1/(e^{βω}-1)`` is
the Bose function.

Intensities can also be calculated for `SampledCorrelations` associated with
classical spin dynamics. In this case, thermal broadening will already be
present, and the line-broadening `kernel` may be omitted. Conversely, the
parameter `kT` becomes required. If positive, it will introduce an intensity
correction factor ``|βω (1 + n_B(ω))|`` that undoes the occupation factor for
the classical Boltzmann distribution and applies the quantum thermal occupation
factor. The special choice `kT = nothing` will suppress the classical-to-quantum
correction factor, and yield statistics consistent with the classical Boltzmann
distribution. If a `kernel` is provided, it will be used to perform a
convolution along the energy axis on top of any intrinsic broadening already
present in the correlation data. In this case, `energies` must be an explicit
list specifying the desired output energies, exactly as for `SpinWaveTheory`.
"""
function intensities(swt::AbstractSpinWaveTheory, qpts; energies, kernel::AbstractBroadening, kT=0, threaded=false)
    return broaden(intensities_bands(swt, qpts; kT, threaded); energies, kernel, threaded)
end

"""
    intensities_static(swt::SpinWaveTheory, qpts; bounds=(-Inf, Inf), kernel=nothing, kT=0, threaded=false)
    intensities_static(sc::SampledCorrelations, qpts; bounds=(-Inf, Inf), kT)
    intensities_static(sc::SampledCorrelationsStatic, qpts)

Like [`intensities`](@ref), but integrates the dynamical correlations
``\\mathcal{S}(𝐪, ω)`` over a range of energies ``ω``. By default, the
integration `bounds` are ``(-∞, ∞)``, yielding the instantaneous (equal-time)
correlations.

In [`SpinWaveTheory`](@ref), the integral will be realized as a sum over
discrete bands. Alternative calculation methods are
[`SpinWaveTheorySpiral`](@ref) and [`SpinWaveTheoryKPM`](@ref).

Classical dynamics data in [`SampledCorrelations`](@ref) can also be used to
calculate static intensities. In this case, the domain of integration will be a
finite grid of available `energies`. Here, the parameter `kT` will be used to
account for the quantum thermal occupation of excitations, as documented in
[`intensities`](@ref).

Static intensities calculated from [`SampledCorrelationsStatic`](@ref) are
dynamics-independent. Instead, instantaneous correlations sampled from the
classical Boltzmann distribution will be reported.
"""
function intensities_static(swt::AbstractSpinWaveTheory, qpts; bounds=(-Inf, Inf), kernel=nothing, kT=0, threaded=false)
    res = intensities_bands(swt, qpts; kT, threaded)  # TODO: with_negative=true
    data_reduced = zeros(eltype(res.data), size(res.data)[2:end])
    for iq in CartesianIndices(data_reduced), ib in axes(res.data, 1)
        ϵ = res.disp[ib, iq]
        if isnothing(kernel) || bounds == (-Inf, Inf)
            if bounds[1] <= ϵ < bounds[2]
                data_reduced[iq] += res.data[ib, iq]
            end
        else
            isnothing(kernel.integral) && error("Kernel must provide integral")
            ihi = kernel.integral(bounds[2] - ϵ)
            ilo = kernel.integral(bounds[1] - ϵ)
            data_reduced[iq] += res.data[ib, iq] * (ihi - ilo)
        end
    end
    StaticIntensities(res.crystal, res.qpts, data_reduced)
end
