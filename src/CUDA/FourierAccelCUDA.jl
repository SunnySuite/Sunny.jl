# TODO
"""
Dipole-dipole interactions computed in Fourier-space. Should produce
identical results (up to numerical precision) as `DipoleReal`, but
is asymptotically faster.
"""
struct DipoleFourierCUDA <: AbstractInteractionCUDA
    int_mat     :: CUDA.CuArray{ComplexF64, 7}
    _spins_ft   :: CUDA.CuArray{ComplexF64, 5}  # Space for Fourier-transforming spins
    _field_ft   :: CUDA.CuArray{ComplexF64, 5}  # Space for holding Fourier-transformed fields
    _field_real :: CUDA.CuArray{Float64, 5}     # Space for holding IFT-transformed fields
    _plan       :: rFTPlan  # Need to become CUFFT types?
    _ift_plan   :: rIFTPlan
end