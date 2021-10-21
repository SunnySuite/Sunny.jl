"Functions for computing energies in Fourier space. All functions expect ±-compressed interaction tensors."

"Computes the Fourier transform of a (spatially compressed) dipole interaction matrix"
function _rfft_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Complex{Float64}}
    A = _reinterpret_dipole_tensor(A)
    rfft(A, 5:ndims(A))
end

"Fourier transforms a dipole system"
function _rfft_dipole_sys(spins::Array{Vec3}) :: Array{Complex{Float64}}
    Sr = reinterpret(reshape, Float64, spins)
    rfft(Sr, 3:ndims(Sr))
end

# FFTW types for various relevant Fourier transform plans using in this file
const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

"""
Dipole-dipole interactions computed in Fourier-space. Should produce
identical results (up to numerical precision) as `DipoleReal`, but
is asymptotically faster.
"""
struct DipoleFourierCPU <: InteractionCPU
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: rFTPlan
    _ift_plan   :: rIFTPlan
end

function DipoleFourierCPU(dip::DipoleDipole, crystal::Crystal, latsize)
    @unpack strength, extent, η = dip
    lattice = Lattice(crystal, latsize)

    A = strength .* precompute_dipole_ewald(lattice; extent=extent, η=η)
    FA = _rfft_dipole_tensor(A)
    nb = nbasis(lattice)
    rftdim = div(size(lattice, 2), 2) + 1
    spins_ft = Array{ComplexF64, 5}(undef, 3, nb, rftdim, size(lattice)[3:end]...)
    field_ft = Array{ComplexF64, 5}(undef, 3, nb, rftdim, size(lattice)[3:end]...)
    field_real = Array{Float64, 5}(undef, 3, size(lattice)...)

    mock_spins = zeros(3, size(lattice)...)
    plan = plan_rfft(mock_spins, 3:ndims(mock_spins); flags=FFTW.MEASURE)
    ift_plan = plan_irfft(spins_ft, size(lattice, 2), 3:ndims(mock_spins); flags=FFTW.MEASURE)
    DipoleFourierCPU(FA, spins_ft, field_ft, field_real, plan, ift_plan)
end

function energy(spins::Array{Vec3, 4}, dip::DipoleFourierCPU)
    FA = dip.int_mat
    FS = dip._spins_ft
    nb = size(spins, 1)
    latsize = size(spins)[2:end]
    even_rft_size = latsize[1] % 2 == 0
    Fsize = size(FS)[3:end]
    spins = _reinterpret_from_spin_array(spins)

    U = 0.0
    mul!(FS, dip._plan, spins)
    # Need to add a normalization factor to some components due to rfft usage
    if even_rft_size
        @views FS[:, :, 2:end-1, :, :] .*= √2
    else
        @views FS[:, :, 2:end, :, :] .*= √2
    end
    @tullio U += real(
        conj(FS[is, ib, j, k, l]) * FA[is, js, ib, jb, j, k, l] * FS[js, jb, j, k, l]
    )
    return U / prod(latsize)
end

"Accumulates the local field coming from dipole interactions, using Fourier transforms"
function _accum_field!(H::Array{Vec3, 4}, spins::Array{Vec3, 4}, dip::DipoleFourierCPU)
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
