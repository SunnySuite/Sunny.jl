"Functions for computing energies in Fourier space. All functions expect ±-compressed interaction tensors."

"Computes the Fourier transform of a (spatially compressed) dipole interaction matrix"
function _rfft_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Complex{Float64}}
    A = _reinterpret_dipole_tensor(A)
    rfft(A, 5:ndims(A))
end

"Fourier transforms a dipole system"
function _rfft_dipole_sys(sys::SpinSystem{D}) :: Array{Complex{Float64}} where {D}
    Sr = reinterpret(reshape, Float64, sys.sites)
    rfft(Sr, 3:ndims(Sr))
end

function DipoleFourier(strength::Float64, lattice::Lattice{3}; extent::Int=4, η::Float64=0.5)
    A = strength .* precompute_dipole_ewald_c(lattice; extent=extent, η=η)
    FA = _rfft_dipole_tensor(A)
    nb = nbasis(lattice)
    rftdim = div(size(lattice, 2), 2) + 1
    spins_ft = Array{ComplexF64, 5}(undef, 3, nb, rftdim, size(lattice)[3:end]...)
    field_ft = Array{ComplexF64, 5}(undef, 3, nb, rftdim, size(lattice)[3:end]...)
    field_real = Array{Float64, 5}(undef, 3, size(lattice)...)

    mock_spins = zeros(3, size(lattice)...)
    plan = plan_rfft(mock_spins, 3:ndims(mock_spins); flags=FFTW.MEASURE)
    ift_plan = plan_irfft(spins_ft, size(lattice, 2), 3:ndims(mock_spins); flags=FFTW.MEASURE)
    DipoleFourier(FA, spins_ft, field_ft, field_real, plan, ift_plan)
end

function energy(sys::SpinSystem{3}, dip::DipoleFourier)
    FA = dip.int_mat
    FS = dip._spins_ft
    nb = nbasis(sys.lattice)
    Fsize = size(FS)[3:end]
    spins = _reinterpret_from_spin_array(sys.sites)
    even_size = sys.lattice.size[1] % 2 == 0

    U = 0.0
    mul!(FS, dip._plan, spins)
    # Need to add a normalization factor to some components due to rfft usage
    if even_size
        @views FS[:, :, 2:end-1, :, :] .*= √2
    else
        @views FS[:, :, 2:end, :, :] .*= √2
    end
    @tullio U += real(
        conj(FS[is, ib, j, k, l]) * FA[is, js, ib, jb, j, k, l] * FS[js, jb, j, k, l]
    )
    return U / prod(sys.lattice.size)
end

"Accumulates the local field coming from dipole interactions, using Fourier transforms"
function _accum_field!(H::Array{Vec3}, spins::Array{Vec3}, dip::DipoleFourier)
    FA = dip.int_mat
    FS = dip._spins_ft
    Fϕ = dip._field_ft
    ϕ  = dip._field_real
    nb = size(spins, 1)
    Fsize = size(FS)[3:end]

    fill!(Fϕ, 0.0)
    spins = _reinterpret_from_spin_array(spins)
    mul!(FS, dip._plan, spins)
    @tullio grad=false Fϕ[s1,b1,i,j,k] += FA[s1,s2,b1,b2,i,j,k] * FS[s2,b2,i,j,k]
    mul!(ϕ, dip._ift_plan, Fϕ)
    ϕ = _reinterpret_to_spin_array(ϕ)
    for i in eachindex(H)
        H[i] = H[i] - 2 * ϕ[i]
    end
end

"Tests these field-using functions give the same answer as `ewald_sum_dipole`"
function test_energy_consistency(sys::SpinSystem{3})
    dip_real = DipoleReal(1.0, sys; extent=4, η=0.5)
    dip_fourier = DipoleFourier(1.0, sys; extent=4, η=0.5)

    direct_energy = ewald_sum_dipole(sys; extent=4, η=0.5)
    real_energy = energy(sys, dip_real)
    fourier_energy = energy(sys, dip_fourier)

    @assert direct_energy ≈ real_energy      "`DipoleRealPre` energy not correct!"
    @assert direct_energy ≈ fourier_energy   "`DipoleFourier` energy not correct!"
end

function test_field_consistency(sys::SpinSystem{3})
    dip_real = DipoleReal(1.0, sys; extent=4, η=0.5)
    dip_fourier = DipoleFourier(1.0, sys; extent=4, η=0.5)

    H1 = zero(sys)
    H2 = zero(sys)
    _accum_field!(H1, sys, dip_real)
    _accum_field!(H2, sys, dip_fourier)

    @assert all(H1 .≈ H2)
end