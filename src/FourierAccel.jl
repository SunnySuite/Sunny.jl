"Functions for computing energies in Fourier space. All functions expect ±-compressed interaction tensors."

"Computes the Fourier transform of a (spatially compressed) dipole interaction matrix"
function _rfft_dipole_tensor(A::OffsetArray{Mat3}) :: Array{Complex{Float64}}
    A = _reinterpret_dipole_tensor(A)
    FFTW.rfft(A, 3:5)
end

"Fourier transforms a dipole system"
function _rfft_dipole_sys(dipoles::Array{Vec3}) :: Array{Complex{Float64}}
    Sr = reinterpret(reshape, Float64, dipoles)
    FFTW.rfft(Sr, 2:4)
end

# FFTW types for various relevant Fourier transform plans using in this file
const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

"""
Dipole-dipole interactions computed in Fourier-space. Should produce
identical results (up to numerical precision) as `DipoleReal`, but
is asymptotically faster.
"""
struct DipoleFourierCPU <: AbstractInteractionCPU
    int_mat     :: Array{ComplexF64, 7}
    _spins_ft   :: Array{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: Array{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: Array{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: rFTPlan
    _ift_plan   :: rIFTPlan
end

function DipoleFourierCPU(dip::DipoleDipole, crystal::Crystal, latsize, site_infos::Vector{SiteInfo};
                          μB=BOHR_MAGNETON::Float64, μ0=VACUUM_PERM::Float64)
    (; extent, η) = dip
    lattice = Lattice(crystal, latsize)

    A = (μ0/4π) * μB^2 .* precompute_dipole_ewald(lattice; extent, η)
    # Conjugate each matrix by the correct g matrices
    for b2 in 1:nbasis(crystal)
        g2 = site_infos[b2].g
        for b1 in 1:nbasis(crystal)
            g1 = site_infos[b1].g
            for ijk in CartesianIndices(axes(A)[1:end-2])
                A[ijk, b1, b2] = g1' * A[ijk, b1, b2] * g2
            end
        end
    end
    FA = _rfft_dipole_tensor(A)
    nb = nbasis(lattice)

    rftdim = div(size(lattice, 1), 2) + 1   
    spins_ft = Array{ComplexF64, 5}(undef, 3, rftdim, size(lattice)[2:3]..., nb)  
    field_ft = Array{ComplexF64, 5}(undef, 3, rftdim, size(lattice)[2:3]..., nb)
    field_real = Array{Float64, 5}(undef, 3, size(lattice)...)

    mock_spins = zeros(3, size(lattice)...)
    plan = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(spins_ft, size(lattice, 1), 2:4; flags=FFTW.MEASURE)

    DipoleFourierCPU(FA, spins_ft, field_ft, field_real, plan, ift_plan)
end

function energy(dipoles::Array{Vec3, 4}, dip::DipoleFourierCPU)
    FA = dip.int_mat
    FS = dip._spins_ft
    latsize = size(dipoles)[1:3]
    even_rft_size = latsize[1] % 2 == 0
    spins = _reinterpret_from_spin_array(dipoles)

    U = 0.0
    mul!(FS, dip._plan, spins)
    # Need to add a normalization factor to some components due to rfft usage
    if even_rft_size
        @views FS[:, 2:end-1, :, :, :] .*= √2
    else
        @views FS[:, 2:end, :, :, :] .*= √2
    end
    @tullio U += real(
        conj(FS[is, j, k, l, ib]) * FA[is, js, j, k, l, ib, jb] * FS[js, j, k, l, jb]
    )
    return U / prod(latsize)
end

"Accumulates the local -∇ℋ coming from dipole-dipole couplings into `B`"
function _accum_neggrad!(H::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, dipdip::DipoleFourierCPU)
    FA = dipdip.int_mat
    FS = dipdip._spins_ft
    Fϕ = dipdip._field_ft
    ϕ  = dipdip._field_real

    fill!(Fϕ, 0.0)
    spins = _reinterpret_from_spin_array(dipoles)
    mul!(FS, dipdip._plan, spins)
    @tullio grad=false Fϕ[s1,i,j,k,b1] += FA[s1,s2,i,j,k,b1,b2] * FS[s2,i,j,k,b2]
    mul!(ϕ, dipdip._ift_plan, Fϕ)
    ϕ = _reinterpret_to_spin_array(ϕ)
    for i in eachindex(H)
        H[i] = H[i] - 2 * ϕ[i]
    end
end
