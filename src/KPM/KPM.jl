# This code will implement the Kernal Polynomial method to approximate the
# dynamical spin structure factor (DSSF). The method takes advantage of the fact
# that the DSSF can be written as a matrix function and hence we can avoid
# diagonalizing the Hamiltonian. Futher, implementing the KPM we do not even
# need to explicity construct the matrix function, instead we only need to
# evaluate the application of the function to some vector. For sparce matrices
# this is efficient making this approximate approach useful in systems with
# large unit cells. 

# Smoothly approximate a Heaviside step function
function regularization_function(ω,σ)
    if ω < 0 
        return 0.0
    elseif 0 ≤ ω ≤ σ
        return (4 - (3ω ./σ)) .*(ω.^3/σ.^3)
    else
        return 1.0
    end
end


"""
    get_all_coefficients(M, ωs, broadening, σ, kT, γ)  

Retrieves the Chebyshev coefficients up to index, M for a user-defined
lineshape. A numerical regularization is applied to treat the divergence of the
dynamical correlation function at small energy. Here, σ² represents an energy
cutoff scale which should be related to the energy resolution. γ is the maximum
eigenvalue used to rescale the spectrum to lie on the interval [-1,1].
Regularization is treated using a cubic cutoff function and the negative
eigenvalues are zeroed out.
"""
function get_all_coefficients(M, ωs, broadening, σ, kT, γ; η=1.0)
    f(ϵ, ω) = regularization_function(ϵ, η*σ) * broadening(ϵ, ω) * thermal_prefactor(kT, ϵ)

    output = zeros(M, length(ωs))
    for i in eachindex(ωs)
        output[:, i] = cheb_coefs(M, 2M, x -> f(γ*x, ωs[i]), (-1, 1))
    end

    return output
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

function kpm_dssf(swt::SpinWaveTheory, qpts, energies, P::Int64, kT, σ, kernel)
    qpts = convert(AbstractQPoints, qpts)
    qs = qpts.qs

    # P is the max Chebyshyev coefficient
    (; sys, measure) = swt
    cryst = orig_crystal(sys)
    qs = Vec3.(qs)
    Nf = nflavors(swt)
    Na = length(eachsite(sys))
    L = Nf*Na
    n_iters = 50
    Avec_pref = zeros(ComplexF64, Na) # initialize array of some prefactors

    Nobs = size(measure.observables, 1)
    Ncorr = length(measure.corr_pairs)
    moments = zeros(ComplexF64, Ncorr, P)
    corrbuf = zeros(ComplexF64, Ncorr)
    intensity = zeros(eltype(measure), length(energies), length(qs))

    u = zeros(ComplexF64, Nobs, 2L) # TODO: (Nobs, Nf, Na, 2)
    α0 = zeros(ComplexF64, Nobs, 2L)
    α1 = zeros(ComplexF64, Nobs, 2L)
    α2 = zeros(ComplexF64, Nobs, 2L)

    for qidx in CartesianIndices(qs)
        q = qs[qidx]
        q_reshaped = to_reshaped_rlu(sys, q)
        q_global = cryst.recipvecs * q

        # Bound eigenvalue magnitude
        lo, hi = eigbounds(swt, q_reshaped, n_iters; extend=0.25)
        γ = max(abs(lo), abs(hi))

        # u(q) calculation)
        for i in 1:Na
            r = sys.crystal.positions[i]
            phase = exp(2π*im * dot(q_reshaped, r))
            Avec_pref[i] = phase / √Na
        end

        if sys.mode == :SUN
            data = swt.data::SWTDataSUN
            N = sys.Ns[1]
            for i in 1:Na, μ in 1:Nobs
                @views O = data.observables_localized[μ, i]
                for f in 1:Nf
                    u[μ, f + (i-1)*Nf]     = Avec_pref[i] * O[f, N] # u[μ, f, i, 1]
                    u[μ, f + (i-1)*Nf + L] = Avec_pref[i] * O[N, f] # u[μ, f, i, 2]
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            data = swt.data::SWTDataDipole
            for i in 1:Na
                sqrt_halfS = data.sqrtS[i]/sqrt(2)
                for μ in 1:Nobs
                    O = data.observables_localized[μ, i]
                    u[μ, i]   = Avec_pref[i] * sqrt_halfS * (O[1] + im*O[2])
                    u[μ, i+L] = Avec_pref[i] * sqrt_halfS * (O[1] - im*O[2])
                end
            end
        end
        
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

        ωdep = get_all_coefficients(P, energies, kernel, σ, kT, γ)

        for iω in eachindex(energies)
            for i in 1:Ncorr
                corrbuf[i] = dot(view(ωdep, :, iω), view(moments, i, :))
            end
            intensity[iω, qidx] = measure.combiner(q_global, corrbuf)
        end
    end

    return BroadenedIntensities(cryst, qpts, energies, intensity)
end

struct SpinWaveTheoryKPM
    swt :: SpinWaveTheory
    P :: Int
    σ :: Float64

    function SpinWaveTheoryKPM(sys::System; P::Int, σ::Float64, measure::Union{Nothing, MeasureSpec}, regularization=1e-8)
        return new(SpinWaveTheory(sys; measure, regularization), P, σ)
    end
end


function intensities!(data, swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, formfactors=nothing, kT=0.0)
    (; swt, P, σ) = swt_kpm
    qpts = convert(AbstractQPoints, qpts)
    @assert size(data) == (length(energies), length(qpts.qs))
    @assert eltype(data) == eltype(swt.measure)

    kpm_intensities(swt, qpts.qs, energies, P, kT, σ, kernel)
    
    return BroadenedIntensities(bands.crystal, bands.qpts, collect(energies), data)
end

function intensities(swt_kpm::SpinWaveTheoryKPM, qpts; energies, kernel::AbstractBroadening, formfactors=nothing, kT=0.0)
    qpts = convert(AbstractQPoints, qpts)
    data = zeros(eltype(swt_kpm.swt.measure), length(energies), length(qpts.qs))
    return intensities!(data, swt_kpm, qpts; energies, kernel, formfactors, kT)
end
