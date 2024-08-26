"""
    SpinWaveTheoryKPM(sys::System; measure, resolution, regularization=1e-8)

**Experimental**

An alternative to [`SpinWaveTheory`](@ref) that uses the kernel polynomial
method (KPM) to perform [`intensities`](@ref) calculations. In traditional spin
wave theory calculations, one would explicitly diagonalize the dynamical matrix,
with a cost that scales like ``𝒪(V^3)`` in the volume ``V`` of the magnetic
cell. KPM instead approximates intensities using polynomial expansion of the
dynamical matrix. The computational cost of KPM scales like ``𝒪(V P)`` in the
polynomial order `P`, and is favorable to direct diagonalization for
sufficiently large magnetic cells.

The polynomial order `P` scales like the spectral bandwidth of the dynamical
matrix divided by the target energy `resolution`. If the specified resolution is
too small (relative to the line broadening kernel), the calculated intensities
will exhibit artificial oscillations in energy.
"""
struct SpinWaveTheoryKPM
    swt :: SpinWaveTheory
    resolution :: Float64
    screening_factor :: Float64

    function SpinWaveTheoryKPM(sys::System; measure::Union{Nothing, MeasureSpec}, regularization=1e-8, resolution, screening_factor=1.0)
        return new(SpinWaveTheory(sys; measure, regularization), resolution, screening_factor)
    end
end


# Smoothly approximate a Heaviside step function
function regularization_function(y)
    if y < 0
        return 0.0
    elseif 0 ≤ y ≤ 1
        return (4 - 3y) * y^3
    else
        return 1.0
    end
end


function mul_Ĩ!(y, x)
    L = size(y, 2) ÷ 2
    view(y, :, 1:L)    .= .+view(x, :, 1:L)
    view(y, :, L+1:2L) .= .-view(x, :, L+1:2L)
end

function mul_A!(swt, y, x, qs_reshaped, γ)
    L = size(y, 2) ÷ 2
    mul_dynamical_matrix!(swt, y, x, qs_reshaped)
    view(y, :, 1:L)    .*= +1/γ
    view(y, :, L+1:2L) .*= -1/γ
end

function set_moments!(moments, measure, u, α)
    map!(moments, measure.corr_pairs) do (μ, ν)
        dot(view(u, μ, :), view(α, ν, :))
    end
end


function intensities!(data, swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, formfactors=nothing, kT=0.0)
    qpts = convert(AbstractQPoints, qpts)

    (; swt, resolution, screening_factor) = swt_kpm
    (; sys, measure) = swt
    cryst = orig_crystal(sys)

    @assert eltype(data) == eltype(measure)
    @assert size(data) == (length(energies), length(qpts.qs))

    Na = length(eachsite(sys))
    Ncells = Na / natoms(cryst)
    Nf = nflavors(swt)
    L = Nf*Na
    n_iters = 50
    Avec_pref = zeros(ComplexF64, Na) # initialize array of some prefactors

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)
    moments = ElasticArray{ComplexF64}(undef, Ncorr, 0)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sys.crystal)

    u = zeros(ComplexF64, Nobs, 2L)
    α0 = zeros(ComplexF64, Nobs, 2L)
    α1 = zeros(ComplexF64, Nobs, 2L)
    α2 = zeros(ComplexF64, Nobs, 2L)

    for (iq, q) in enumerate(qpts.qs)
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        # Represent each local observable A(q) as a complex vector u(q) that
        # denotes a linear combination of HP bosons.

        for i in 1:Na
            r = sys.crystal.positions[i]
            Avec_pref[i] = exp(2π*im * dot(q_reshaped, r))
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], norm2(q_global))
        end

        if sys.mode == :SUN
            data_sun = swt.data::SWTDataSUN
            N = sys.Ns[1]
            for i in 1:Na, μ in 1:Nobs
                @views O = data_sun.observables_localized[μ, i]
                for f in 1:Nf
                    u[μ, f + (i-1)*Nf]     = Avec_pref[i] * O[f, N]
                    u[μ, f + (i-1)*Nf + L] = Avec_pref[i] * O[N, f]
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_large_s)
            data_dip = swt.data::SWTDataDipole
            for i in 1:Na
                sqrt_halfS = data_dip.sqrtS[i]/sqrt(2)
                for μ in 1:Nobs
                    O = data_dip.observables_localized[μ, i]
                    u[μ, i]   = Avec_pref[i] * sqrt_halfS * (O[1] + im*O[2])
                    u[μ, i+L] = Avec_pref[i] * sqrt_halfS * (O[1] - im*O[2])
                end
            end
        end

        # Bound eigenvalue magnitudes and determine order of polynomial
        # expansion

        lo, hi = eigbounds(swt, q_reshaped, n_iters; extend=0.25)
        γ = max(abs(lo), abs(hi))
        P = max(round(Int, π*γ/2resolution), 2)
        resize!(moments, Ncorr, P)
        σ = resolution * screening_factor

        # Perform Chebyshev recursion

        q_repeated = fill(q_reshaped, Nobs)
        mul_Ĩ!(α0, u)
        mul_A!(swt, α1, α0, q_repeated, γ)
        set_moments!(view(moments, :, 1), measure, u, α0)
        set_moments!(view(moments, :, 2), measure, u, α1)
        for m in 3:P
            mul_A!(swt, α2, α1, q_repeated, γ)
            @. α2 = 2*α2 - α0
            set_moments!(view(moments, :, m), measure, u, α2)
            (α0, α1, α2) = (α1, α2, α0)
        end

        # Transform Chebyshev moments to intensities for each ω

        buf = zeros(2P)
        plan = FFTW.plan_r2r!(buf, FFTW.REDFT10)

        for (iω, ω) in enumerate(energies)
            f(x) = regularization_function(x / σ) * kernel(x, ω) * thermal_prefactor(x; kT)

            coefs = cheb_coefs!(P, f, (-γ, γ); buf, plan)
            for i in 1:Ncorr
                corrbuf[i] = dot(coefs, view(moments, i, :)) / Ncells
            end
            data[iω, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return Intensities(cryst, qpts, collect(energies), data)
end

function intensities(swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, formfactors=nothing, kT=0.0)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_kpm.swt.measure), length(energies), length(qpts.qs))
    return intensities!(data, swt_kpm, qpts; energies, kernel, formfactors, kT)
end
