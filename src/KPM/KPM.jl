# This code will implement the Kernal Polynomial method to approximate the
# dynamical spin structure factor (DSSF). The method takes advantage of the fact
# that the DSSF can be written as a matrix function and hence we can avoid
# diagonalizing the Hamiltonian. Futher, implementing the KPM we do not even
# need to explicity construct the matrix function, instead we only need to
# evaluate the application of the function to some vector. For sparce matrices
# this is efficient making this approximate approach useful in systems with
# large unit cells. 

"""
    bose_function(kT, x)

Returns the Bose occupation factor for energy, x, at temperature kT.
"""
function bose_function(kT, x) 
    r =  1 / (exp(x/kT) - 1)
    if isinf(r)
      return 0.
    end
    return r
end

"""
    regularization_function(ω,σ)

Returns a regularization factor to apply to the intensity at low energy according to a smooth approximation to a step function
with width, σ. 
"""
function regularization_function(ω,σ)
    if ω < 0 
        return 0.0
    elseif 0 ≤ ω ≤ σ
        return (4 - (3ω ./σ)) .*(ω.^3/σ.^3)
    else
        return 1.0
    end
end

function srf(ω,σ)
    if -σ ≤ ω ≤ σ
        return (3(ω ./σ).^4 .- 2(ω ./σ).^6)
    else
        return 1.0
    end
end



"""
    get_all_coefficients(M, ωs, broadening, σ, kT,γ)  

Retrieves the Chebyshev coefficients up to index, M for a user-defined
lineshape. A numerical regularization is applied to treat the divergence of the
dynamical correlation function at small energy. Here, σ² represents an energy
cutoff scale which should be related to the energy resolution. γ is the maximum
eigenvalue used to rescale the spectrum to lie on the interval [-1,1].
Regularization is treated using a cubic cutoff function and the negative
eigenvalues are zeroed out.
"""
function get_all_coefficients(M, ωs, broadening, σ, kT,γ;η=1.0, regularization_style)
    f = if regularization_style == :cubic
      (ω,x) -> regularization_function(γ*x,η*σ) * broadening(ω, x*γ, σ) * (1 + bose_function(kT, x*γ))
    elseif regularization_style == :tanh
      (ω,x) -> tanh((γ*x/(η*σ))^2) * broadening(ω, x*γ, σ) * (1 + bose_function(kT, x*γ))
    elseif regularization_style == :none
      (ω,x) -> broadening(ω, x*γ, σ) * (1 + bose_function(kT, x*γ))
    elseif regularization_style == :susceptibility
      (ω,x) -> broadening(ω, x*γ, σ)
    elseif regularization_style == :srf
      (ω,x) -> srf(γ*x,η*σ) * broadening(ω, x*γ, σ) * (1 + bose_function(kT, x*γ))
    end
    output = OffsetArray(zeros(M, length(ωs)), 0:M-1, 1:length(ωs))
    for i in eachindex(ωs)
        output[:, i] = cheb_coefs(M, 2M, x -> f(ωs[i], x), (-1, 1))
    end
    return output
end


"""
    get_all_coefficients_legacy(M, ωs, broadening, σ, kT,γ)  

Retrieves the Chebyshev coefficients up to index, M for a user-defined
lineshape. A numerical regularization is applied to treat the divergence of the
dynamical correlation function at small energy. Here, σ² represents an energy
cutoff scale which should be related to the energy resolution. γ is the maximum
eigenvalue used to rescale the spectrum to lie on the interval [-1,1].
Regularization is treated using a tanh cutoff function and the negative
eigenvalues are not zeroed out.
"""
function get_all_coefficients_legacy(M, ωs, broadening, σ, kT,γ)
    f(ω, x) = tanh((x/σ)^2) * broadening(ω, x*γ, σ) * (1 + bose_function(kT, x))
    output = OffsetArray(zeros(M, length(ωs)), 0:M-1, 1:length(ωs))
    for i in eachindex(ωs)
        output[:, i] = cheb_coefs(M, 2M, x -> f(ωs[i], x), (-1, 1))
    end
    return output
end

"""
    apply_kernel(offset_array, kernel, M)

Applies the Jackson kernel (defined in Chebyshev.jl) to an OffsetArray. The
Jackson kernel damps the coefficients of the Chebyshev expansion to reduce
"Gibbs oscillations" or "ringing" -- numerical artifacts present due to the
truncation offset_array the Chebyshev series. It should be noted that employing
the Jackson kernel comes at a cost to the energy resolution.
"""

#  
function apply_kernel(offset_array, kernel, M)
    kernel === "jackson" ? offset_array .= offset_array .* jackson_kernel(M) : 
    return offset_array
end

"""
    kpm_dssf(swt::SpinWaveTheory, qs,ωlist,P::Int64,kT,σ,broadening; kernel)

Calculated the Dynamical Spin Structure Factor (DSSF) using the Kernel Polynomial Method (KPM). Requires input in the form of a 
SpinWaveTheory which contains System information and the rotated spin matrices. The calculation is carried out at each wavevectors
in qs for all energies appearing in ωlist. The Chebyshev expansion is taken to P terms and the lineshape is specified by the user-
defined function, broadening. kT is required for the calcualation of the bose function and σ is the width of the lineshape and 
defines the low energy cutoff σ². There is a keyword argument, kernel, which speficies a damping kernel. 
"""
 
function kpm_dssf(swt::SpinWaveTheory, qs, ωlist, P::Int64, kT, σ, broadening; kernel, regularization_style)
    # P is the max Chebyshyev coefficient
    (; sys, data) = swt
    qs = Vec3.(qs)
    Nm, Ns = length(sys.dipoles), sys.Ns[1] # number of magnetic atoms and dimension of Hilbert space
    Nf = sys.mode == :SUN ? Ns-1 : 1
    N=Nf+1
    nmodes = Nf*Nm
    M = sys.mode == :SUN ? 1 : (Ns-1) # scaling factor (=1) if in the fundamental representation
    sqrt_M = √M #define prefactors
    sqrt_Nm_inv = 1.0 / √Nm #define prefactors
    S = (Ns-1) / 2
    sqrt_halfS  = √(S/2) #define prefactors   
    Ĩ = Diagonal([ones(nmodes); -ones(nmodes)]) 
    n_iters = 50
    Avec_pref = zeros(ComplexF64, Nm) # initialize array of some prefactors   
    chebyshev_moments = OffsetArray(zeros(ComplexF64, 3, 3, length(qs), P), 1:3, 1:3, 1:length(qs), 0:P-1)
    Sαβs = zeros(ComplexF64,3,3,length(qs),length(ωlist))
    for qidx in CartesianIndices(qs)
        q = qs[qidx]
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        u = zeros(ComplexF64,3,2*nmodes)
        lo, hi = eigbounds(swt, q_reshaped, n_iters; extend=0.25) # calculate bounds
         # Upper bound for generalized eigenvalues. Factor of 2 accounts for
         # implicit rescaling of Hamiltonian.
        γ = max(abs(lo), abs(hi))

        # u(q) calculation)
        for site in 1:Nm
            # note that d is the chemical coordinates
            chemical_coor = sys.crystal.positions[site] # find chemical coords
            phase = exp(2*im * π  * dot(q_reshaped, chemical_coor)) # calculate phase
            Avec_pref[site] = sqrt_Nm_inv * phase  # define the prefactor of the tS matrices
        end

        # calculate u(q)
        if sys.mode == :SUN
            for site=1:Nm
                @views tS_μ = data.observable_operators[:, :, :, site]*Avec_pref[site] 
                for μ=1:3
                    for j=2:N
                        u[μ,(j-1)+(site-1)*(N-1)] = tS_μ[j,1,μ] 
                        u[μ,(N-1)*Nm+(j-1)+(site-1)*(N-1)] = tS_μ[1,j,μ]
                    end
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            for site in 1:Nm
                R = data.local_rotations[site]
                u[1,site]= Avec_pref[site] * sqrt_halfS * (R[1,1] + 1im * R[1,2])  
                u[1,site+nmodes] = Avec_pref[site] * sqrt_halfS * (R[1,1] - 1im * R[1,2])
                u[2,site]= Avec_pref[site] * sqrt_halfS * (R[2,1] + 1im * R[2,2]) 
                u[2,site+nmodes] = Avec_pref[site] * sqrt_halfS * (R[2,1] - 1im * R[2,2]) 
                u[3,site]= Avec_pref[site] * sqrt_halfS * (R[3,1] + 1im * R[3,2]) 
                u[3,site+nmodes] = Avec_pref[site] * sqrt_halfS * (R[3,1] - 1im * R[3,2]) 
            end
        end

        for β in 1:3
            α0 = zeros(ComplexF64,2*nmodes)
            α1 = zeros(ComplexF64,2*nmodes)
            mul!(α0, Ĩ, u[β,:]) # calculate α0
            multiply_by_hamiltonian!(swt, α1, α0, q_reshaped)
            mul!(α1, Ĩ, α1/γ)
            for α in 1:3
                chebyshev_moments[α,β,qidx,0] = dot(u[α,:], α0) #removed symmetrization
                chebyshev_moments[α,β,qidx,1] = dot(u[α,:], α1) #removed symmetrization
            end
            for m in 2:P-1
                αnew = zeros(ComplexF64, 2*nmodes)
                multiply_by_hamiltonian!(swt, αnew, α1, q_reshaped)
                mul!(αnew, Ĩ, αnew/γ)
                @. αnew = 2*αnew - α0
                for α in 1:3
                    chebyshev_moments[α, β, qidx, m] = dot(u[α,:], αnew) #removed symmetrization
                end
                (α1, α0) = (αnew, α1)
            end
        end
        ωdep = get_all_coefficients(P, ωlist, broadening, σ, kT, γ; regularization_style)
        apply_kernel(ωdep, kernel, P)
        for w in eachindex(ωlist)
            for α in 1:3, β in 1:3
                Sαβs[α, β, qidx, w] = sum(chebyshev_moments[α, β, qidx,:] .*  ωdep[:,w])
            end
        end 
    end
    return Sαβs
end

"""
    kpm_intensities(swt::SpinWaveTheory, qs,ωlist,P::Int64,kT,σ,broadening; kernel)

Calculated the neutron scattering intensity using the Kernel Polynomial Method (KPM). Calls KPMddsf and so takes the same parameters.
Requires input in the form of a SpinWaveTheory which contains System information and the rotated spin matrices. The calculation is 
carried out at each wavevectors in qs for all energies appearing in ωlist. The Chebyshev expansion is taken to P terms and the 
lineshape is specified by the user-defined function, broadening. kT is required for the calcualation of the bose function and σ is 
the width of the lineshape and defines the low energy cutoff σ². There is an optional keyword argument, kernel, which speficies a 
damping kernel. The default is to include no damping.  
"""
function kpm_intensities(swt::SpinWaveTheory, qs, ωvals,P::Int64,kT,σ,broadening; kernel = nothing,regularization_style)
    (; sys) = swt
    qs = Vec3.(qs)
    Sαβs = kpm_dssf(swt, qs, ωvals, P, kT, σ, broadening; kernel, regularization_style)
    num_ω = length(ωvals)
    is = zeros(Float64, size(qs)..., num_ω)
    for qidx in CartesianIndices(qs)
        q_reshaped = to_reshaped_rlu(sys, qs[qidx])
        q_absolute = sys.crystal.recipvecs * q_reshaped
        polar_mat = polarization_matrix(q_absolute)
        is[qidx, :] = real(sum(polar_mat .* Sαβs[:,:,qidx,:],dims=(1,2)))
    end
    return is
end


#=
function Itilde!(α, n)
    view(α, n+1:2n) .*= -1
end

function Run_Recurrence_fast(swt,q_reshaped,γ,u,nmodes,chebyshev_moments,M)
    α0 = zeros(ComplexF64,2*nmodes)
    α1 = zeros(ComplexF64,2*nmodes)
    α2 = zeros(ComplexF64,2*nmodes)

    chebyshev_moments .= 0
    for β=1:3
        # α0 = ̃I u
        α0 .= view(u, β, :)
        Itilde!(α0, nmodes)

        # α1 = ̃I A α0
        Sunny.multiply_by_hamiltonian_dipole!(α1,α0,swt,q_reshaped)
        @. α1 = α1 * (2/γ)
        Itilde!(α1, nmodes)

        for α=1:3
            chebyshev_moments[α,β,0] = dot(view(u, α,:), α0)
            chebyshev_moments[α,β,1] = dot(view(u, α,:), α1)
        end
        
        for m=2:M-1
            # α2 = ̃2 I A α1 - α0
            Sunny.multiply_by_hamiltonian_dipole!(α2,α1,swt,q_reshaped) 
            @. α2 = α2 * (2/γ)
            Itilde!(α2, nmodes)
            @. α2 = 2*α2 - α0

            for α=1:3
                chebyshev_moments[α,β,m] = dot(view(u, α,:),α2)
            end
            (α1, α0, α2) = (α2, α1, α0)
        end
    end
end
=#
