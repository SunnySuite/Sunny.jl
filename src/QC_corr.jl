# * Correction for Quantum-Classical Correspondence ver 1
# * Written by Chaebin Kim 09/20/2022
# * Ref) S. Zhang et al. Physical Review Letters 122, 167203 (2019)

# ? QC_corr! : Calculating Quantum-Classical correspondence βω
# ? hexa_corr! : Correcting the momentum space with hexagonal unit cell
# ? Cobalt_ff! : Calculating magnetic form factor for Co2+ correctly
# TODO: Currently, the output of hexa_corr cannot replace the MySQWperp due to the different size of OffsetArray...
# TODO: This could be better than this

"""
    QC_corr!(MySQWperp, EN, kT, Nω)

Multiplying the Quantum-Classical Correspondence factor βω as
'''math
    S(q,ω)_{Quantum} = βωS(q,ω)_{Classical}
'''
where β = 1/kT and ω is energy.

"""
function QC_corr!(MySQWperp, EN, kT, Nω)
    β = 1 / kT
    βω = reshape(β * EN, 1, 1, 1, Nω) # Correction factor βω
    βω = repeat(βω, outer=[size(MySQWperp.sfactor)[1], size(MySQWperp.sfactor)[2], size(MySQWperp.sfactor)[3], 1])
    Ls = MySQWperp.sfactor.offsets .+ 1
    Le = -1 .* MySQWperp.sfactor.offsets
    QC_corr = OffsetArray(βω, Ls[1]:Le[1], Ls[2]:Le[2], Ls[3]:Le[3], Ls[4]:Le[4])

    MySQWperp.sfactor .= MySQWperp.sfactor .* QC_corr
    MySQWperp.sfactor .= MySQWperp.sfactor ./ maximum(MySQWperp.sfactor) # Renormalization!

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