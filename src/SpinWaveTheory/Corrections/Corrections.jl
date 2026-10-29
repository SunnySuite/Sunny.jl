# Corrections to linear spin wave theory (LSWT) at sub-leading order in 1/s.
#
# Expanding the spin Hamiltonian in Holstein-Primakoff bosons b, b† yields
#
#     H = E_cl + H₁ + H₂ + H₃ + H₄ + …,
#
# where Hₙ collects the n-boson terms and scales as s^(2-n/2). At a classical
# energy minimum H₁ vanishes. LSWT retains H₂ only, which `bogoliubov!`
# diagonalizes into quasi-particle modes α, α† with energies ε. Relative to
# LSWT, both H₄ (at first order) and H₃ (at second order) contribute at O(1/s).
#
# Writing quasi-particle momenta as subscripts and letting N denote the number
# of magnetic cells, the cubic term splits into a decay and a source channel,
#
#     H₃ = (1/2!√N) Σ_{1+2=3} [Γ₁(3; 1, 2) α†₃ α₁ α₂ + h.c.]
#        + (1/3!√N) Σ_{1+2+3=0} [Γ₂(1, 2, 3) α†₁ α†₂ α†₃ + h.c.],
#
# whose self-energies SelfEnergy.jl writes down. At T = 0 only the decay channel
# resonates at ω > 0, giving magnons a finite lifetime.
#
# Sunny's Nambu conventions are reused throughout. In particular the columns
# `L+1:2L` of a Bogoliubov matrix `T` obtained at wavevector `q` are the
# eigenvectors at `-q` (see `excitations!`), so a band index ranging over the
# full Nambu space `1:2L` reaches both channels at once.
#
# `bogoliubov!` fixes the phase of each band independently, so nothing ties the
# eigenvectors of one call to those of another, and degenerate bands may even
# get different bases. Every expression that combines amplitudes is therefore
# built from a single diagonalization per wavevector. Where T(-𝐪) is needed
# alongside T(𝐪), it comes from the exact Nambu symmetry P H(𝐪)ᵗ P = H(-𝐪),
# with P the swap a ↦ nambu_conj(a, L), rather than from a second call. The
# symmetry follows from H22(𝐪) = H11(-𝐪)ᵗ and H12(𝐪) = H12(-𝐪)ᵗ.
#
# The expansion is about the Gaussian state of `swt`, its vacuum; see Vacuum.jl.
# Static corrections are contractions in it (StaticCorrections.jl), the cubic
# vertex gives a frequency-dependent self-energy (SelfEnergy.jl), and
# CorrectedIntensities.jl resums both into the structure factor. The test suite
# certifies every term by comparing to exact diagonalization of a cluster with
# anisotropic interactions and readouts.
#
#     Vacuum       to_hartree_fock_corrected_vacuum, to_self_consistent_vacuum,
#                  GoldstonePin
#     Spectra      corrected_intensities
#     Statics      corrected_energy_per_site, gaussian_energy_per_site,
#                  boson_density, corrected_magnetic_moments
#     Pieces       hartree_fock_correction, tadpole_correction,
#                  observable_corrections, anisotropy_correction,
#                  static_self_energy, cubic_self_energy
#     Internals    OneLoop, SelfEnergy, corrected_channels
#
# Their keyword `bz` selects the Brillouin-zone grid, and the static corrections
# take `tol` for adaptive cubature; see `BZIntegration`.

# Errors unless `swt` describes a model for which the corrections in this
# directory are implemented
function check_corrections_supported(swt::SpinWaveTheory)
    isnothing(swt.sys.ewald) || error("Corrected spin wave theory does not yet support Ewald interactions")
end
