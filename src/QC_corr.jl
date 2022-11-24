# * Correction for Quantum-Classical Correspondence
# *
# * Written by Chaebin Kim 09/20/2022
# * revised by Chaebin Kim 11/23/2022
# *
# * Ref) S. Zhang et al. Physical Review Letters 122, 167203 (2019)
# * Ref) E. M. Smith et al. Physical Review X 12, 021015 (2022)

# ? QC_corr! : Calculating Quantum-Classical correspondence
# ? hexa_corr! : Correcting the momentum space with hexagonal unit cell
# ? Cobalt_ff! : Calculating magnetic form factor for Co2+ correctly
# TODO: Currently, the output of hexa_corr cannot replace the MySQWperp due to the different size of OffsetArray...
# TODO: This could be better than this

"""
    QC_corr!(MySQWperp, kT)

Multiplying the Quantum-Classical Correspondence factor βω/(1-e^(-βω)) as
'''math
    S(q,ω)_{Quantum} = βω/(1-e^(-βω))S(q,ω)_{Classical}
'''
where β = 1/kT and ω is energy.

"""
function QC_corr!(MySQWperp, kT)
    N1, N2, N3, N4 = size(MySQWperp.sfactor)
    Ls = MySQWperp.sfactor.offsets .+ 1
    Le = -1 .* MySQWperp.sfactor.offsets

    β = 1 / kT
    wmax = round(2 * pi / (MySQWperp.dt * MySQWperp.meas_period))
    EN = range(Ls[4], length=N4) * 2 * wmax / N4
    
    tmp = β * EN ./ (ones(size(EN)) - exp.(-1 * β * EN))
    tmp[round(Int, N4 / 2)] = 1.0 # remove divergence value at E = 0

    βω = reshape(tmp, 1, 1, 1, N4) # Correction factor βω
    βω = repeat(βω, outer=[N1, N2, N3, 1])

    CQ_corr = OffsetArray(βω, Ls[1]:Le[1], Ls[2]:Le[2], Ls[3]:Le[3], Ls[4]:Le[4])

    MySQWperp.sfactor .= MySQWperp.sfactor .* CQ_corr

    return MySQWperp
end

"""
    hexa_corr!(MySQWperp)

Transformation for hexagonal lattice

"""

function hexa_corr!(MySQWperp)
    A = OffsetArrays.no_offset_view(MySQWperp.sfactor)
    Lx, Ly, Lz, T = size(A)
    epsil = 1e-6 # infinitesimal number for correct round function

    New = zeros(round(Int, 3Lx / 2), Ly, Lz, T)
    for i = 1:Lx
        for j = 1:Ly
            y_pot = round(Int, (i + 2 * j) / 2 + epsil)
            New[y_pot, i, :, :] .= A[j, i, :, :]
        end
    end

    Ls = MySQWperp.sfactor.offsets .+ 1
    Le = -1 .* MySQWperp.sfactor.offsets
    L1s = -1 * round(Int, size(New, 1) / 2) + 1
    L1e = round(Int, size(New, 1) / 2)
    New_OSA = OffsetArray(New, L1s:L1e, Ls[2]:Le[2], Ls[3]:Le[3], Ls[4]:Le[4])

    # MySQWperp.sfactor .= New_OSA

    return New_OSA
end

function Cobalt_ff!(SQW, Q1, Q2, Q3)
    A = OffsetArrays.no_offset_view(SQW)
    Lx, Ly, Lz, T = size(A)

    # Magnetic form factor for cobalt 2+
    fj0 = [0.4332; 14.3553; 0.5857; 4.6077; -0.0382; 0.1338; 0.0179]
    fj2 = [1.9049; 11.6444; 1.3159; 4.3574; 0.3146; 1.6453; 0.0017]

    norm = 1 / (4 ) / pi
    g = 4

    fp = zeros(Lx, Ly, Lz)

    for i = 1:Lx
        for j = 1:Ly
            for k = 1:Lz

                Qq = [Q1[i] * 2 * pi / 3.96; Q2[j] * 2 * pi / 3.96; Q3[k] * 2 * pi / 6.66]
                Qsq = Qq' * Qq
                j0 = fj0[1] * exp(-1 * fj0[2] * Qsq * norm^2) + fj0[3] * exp(-1 * fj0[4] * Qsq * norm^2) + fj0[5] * exp(-1 * fj0[6] * Qsq * norm^2) + fj0[7]
                j2 = (fj2[1] * exp(-1 * fj2[2] * Qsq * norm^2) + fj2[3] * exp(-1 * fj2[4] * Qsq * norm^2) + fj2[5] * exp(-1 * fj2[6] * Qsq * norm^2) + fj2[7]) .* Qsq * norm^2

                fp[i, j, k] = j0 + (2 - g) / g * j2


            end
        end
    end

    fp = reshape(fp, Lx, Ly, Lz, 1)
    fp = repeat(fp, outer=[1, 1, 1, T])

    A .= A .* fp

    Ls = SQW.offsets .+ 1
    Le = -1 .* SQW.offsets
    New_OSA = OffsetArray(A, Ls[1]:Le[1], Ls[2]:Le[2], Ls[3]:Le[3], Ls[4]:Le[4])

    # MySQWperp.sfactor .= New_OSA

    return New_OSA, fp
end