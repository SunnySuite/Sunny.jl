# Boson words of the observables that `measure` reads out, and of the spin
# components. The one-boson words are the LSWT amplitudes; the others feed the
# 1/s corrections of Corrections/ObservableCorrections.jl.

# Monomials of `K` bosons in the expansion of part p of observable μ at site i,
# in the local frame and without the Fourier phase or form factor that
# `observable_prefactor` supplies. The one-boson words are the LSWT amplitudes
# that `observable_coefficients` compiles. Of the others, only K = 2, the
# longitudinal word, and K = 3, the leading correction to the transverse ones,
# are needed.
#
# In mode :SUN an observable is a matrix on the same footing as a term of the
# Hamiltonian, so its words are those of `local_monomials`. In the dipole modes it
# is a spin component, whose words of at most four bosons the truncated
# expansion holds exactly, or a general polynomial in the spin. An observable
# never enters the Hamiltonian, so the truncation rule that protects its
# Goldstone modes does not apply, and the words of a polynomial are those of the
# operator itself, `exact_stevens_words`, at the tree level and in the 1/s
# corrections alike. In mode :dipole_uncorrected the operator is the limit s → ∞
# of its Stevens polynomial, whose words are read off at their leading power,
# and also at the next for the 1/s corrections. Either way the stored
# coefficients are those of the operator, without `rcs_factors`.
function observable_monomials(sys::System, data, μ, i, ::Val{K}; part=1, beyond_tree=true) where K
    L = nbands(sys)
    if sys.mode == :SUN
        A = (data::SWTDataSUN).observables[μ, i, part]
        return local_monomials(A, i, zero(Vec3), Val{K}(), nflavors(sys), L)
    end

    @assert sys.mode in (:dipole, :dipole_uncorrected)
    part == 1 || error("Distinct parts unsupported in dipole-mode.")
    (; sqrtS, observables) = data::SWTDataDipole
    O = observables[μ, i]
    words = if O isa Vec3
        spin_component_words(O, √2 * sqrtS[i])
    elseif sys.mode == :dipole
        exact_stevens_words(stevens_ranks(O), sqrtS[i]^2)
    else
        (lead, sub) = stevens_words(stevens_ranks(O), √2 * sqrtS[i], false)
        beyond_tree ? lead + sub : lead
    end
    return local_word_monomials(words, i, L, Val{K}())
end

# The same for a measure of a single part, which is all that the 1/s corrections
# support
function observable_monomials(swt::SpinWaveTheory, μ, i, ::Val{K}) where K
    @assert num_parts_per_unit(swt.measure) == 1  # Entangled units are rejected
    return observable_monomials(swt.sys, swt.data, μ, i, Val{K}())
end

# Every K-boson word of observable μ, over the sites of the magnetic cell.
observable_monomials(swt::SpinWaveTheory, μ, ::Val{K}) where K =
    reduce(vcat, observable_monomials(swt, μ, i, Val{K}()) for i in 1:nsites(swt.sys))

# The monomials of K bosons of every observable. Those of K = 3 and 4 are what a
# correction to the one-magnon and two-magnon amplitudes contract.
all_observable_monomials(swt::SpinWaveTheory, ::Val{K}) where K =
    reduce(vcat, (observable_monomials(swt, μ, Val{K}()) for μ in 1:num_observables(swt.measure)); init=BosonMonomial{K}[])

# In mode :SUN the spin components are matrices on the same footing as a term of
# the Hamiltonian, so their words are those of `local_monomials`. In the dipole
# modes the Holstein-Primakoff expansion of ExpansionDipole.jl gives the words
# directly, ⟨𝐒⟩ = R (σ Re v, σ Im v, s - ⟨b†b⟩) in a local frame whose ẑ is the
# classical dipole.
function spin_monomials(swt::SpinWaveTheory, α, i, ::Val{K}) where K
    (; sys, data) = swt
    L = nbands(swt)
    o = zero(Vec3)

    if sys.mode == :SUN
        return local_monomials((data::SWTDataSUN).spin_ops[α, i], i, o, Val{K}(), nflavors(swt), L)
    end

    @assert sys.mode in (:dipole, :dipole_uncorrected)
    (; local_rotations, sqrtS) = data::SWTDataDipole
    # Row α of R is the global component α in the local frame. The words are the
    # classical dipole, the tilt, and the depletion. All vanish on a vacant site,
    # where R does.
    words = spin_component_words(local_rotations[i][α, :], √2 * sqrtS[i])
    return local_word_monomials(words, i, L, Val{K}())
end
