"""
    energy_per_site_lswt_correction(swt::SpinWaveTheory; opts...)

Computes the [𝒪(1/λ) or 𝒪(1/S)] correction to the classical energy **per
site** [𝒪(λ²) or 𝒪(S²)] given a [`SpinWaveTheory`](@ref). The correction
[𝒪(λ) or 𝒪(S)] includes a uniform term (For instance, if the classical energy
is αJS², the LSWT gives a correction like αJS) and the summation over the
zero-point energy for all spin-wave modes, i.e., 1/2 ∑ₙ ∫d³q ω(q, n), where q
belongs to the first magnetic Brillouin zone and n is the band index.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the HCubature package documentation
for details.
"""
function energy_per_site_lswt_correction(swt::SpinWaveTheory; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    Natoms = natoms(sys.crystal)
    L = nbands(swt)
    # Create matrix and vector buffers to reuse them
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)
    E = zeros(L)

    if sys.mode == :SUN
        hamiltonian_function! = swt_hamiltonian_SUN!
    else
        @assert sys.mode in (:dipole, :dipole_large_S)
        hamiltonian_function! = swt_hamiltonian_dipole!
    end

    # The uniform term (trace of the (1,1)-block of the spin-wave Hamiltonian)
    q₀ = Vec3(0.0, 0.0, 0.0)
    hamiltonian_function!(H, swt, q₀)
    δE = -real(tr(view(H, 1:L, 1:L))) / Natoms

    # Integrate zero-point energy over the first magnetic Brillouin zone 𝐪 ∈ [0, 1]³ (in RLU)
    res = hcubature((0,0,0), (1,1,1); opts...) do q
        q = Vec3(q)
        # Clear the energy buffer
        E .= 0.0
        hamiltonian_function!(H, swt, q)
        E .= bogoliubov!(V, H)
        return sum(E) / 2Natoms
    end

    println("Zero-point energy is ", res[1], "±", res[2])

    δE += res[1]

    return δE

end

# Calculates the magnetization reduction for :SUN mode for site `i`
function magnetization_lswt_correction_sun(swt::SpinWaveTheory, i::Int; opts...)
    (; sys, data) = swt

    N = sys.Ns[1]
    Natoms = natoms(sys.crystal)
    L = (N - 1) * Natoms

    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)
    S_buf = zeros(ComplexF64, N, N)

    @views dipole = sys.dipoles[i]
    n = normalize(dipole)
    for μ in 1:3
        @views O = data.observable_operators[:, :, μ, i]
        @. S_buf += n[μ] * O
    end

    res = hcubature((0,0,0), (1,1,1); opts...) do q
        q = Vec3(q)
        swt_hamiltonian_SUN!(H, swt, q)
        bogoliubov!(V, H)
        acc = 0.0
        for band in L+1:2L
            v = reshape(view(V, :, band), N-1, Natoms, 2)
            for α in 1:N-1
                for β in 1:N-1
                    acc += -(S_buf[N, N]*δ(α, β) - S_buf[α, β]) * conj(v[α, i, 1]) * v[β, i, 1]
                end
            end
        end
        return real(acc)
    end

    println("Site ", i)
    println("Classical magnetization")
    println(real(S_buf[N, N]))
    println("Correction from LSWT ")
    println(res[1], "±", res[2])
    return res[1]
end

# Calculates the magnetization reduction for :dipole mode for site `i`
function magnetization_lswt_correction_dipole(swt::SpinWaveTheory, i::Int; opts...)
    (; sys) = swt
    N = sys.Ns[1]
    S = (N-1)/2

    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    res = hcubature((0,0,0), (1,1,1); opts...) do q
        swt_hamiltonian_dipole!(H, swt, Vec3(q))
        bogoliubov!(V, H)
        return -norm2(view(V, L+i, 1:L))
    end

    println("Site ", i)
    println("Classical magnetization")
    println(S)
    println("Correction from LSWT ")
    println(res[1], "±", res[2])

    return res[1]
end

"""
    magnetization_lswt_correction(swt::SpinWaveTheory, i::Int; opts...)

Calculates the reduction in the classical magnetization given a
[`SpinWaveTheory`](@ref) from LSWT for magnetic sublattice `i`. In the case of
`:dipole` and `:dipole_large_S` mode, the classical magnetization is always
maximized to spin size `S`. While in `:SUN` mode, the classical magnetization
can be smaller than `S` due to anisotropic interactions.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the HCubature package documentation
for details.
"""
function magnetization_lswt_correction(swt::SpinWaveTheory, i::Int; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    if sys.mode == :SUN
        δS = magnetization_lswt_correction_sun(swt, i; opts...)
    else
        @assert sys.mode in (:dipole, :dipole_large_S)
        δS = magnetization_lswt_correction_dipole(swt, i; opts...)
    end
    return δS
end
