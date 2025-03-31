# These functions mirror closely those of regular SpinWaveTheory. The main
# difference is careful selection of summation over branches (1, 2, 3)
# corresponding to momentum transfers (q-k, q, q+k).

function energy_per_site_lswt_correction(sswt::SpinWaveTheorySpiral; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; swt) = sswt
    Natoms = natoms(swt.sys.crystal)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # Uniform correction at q = 0, branch 2 is the "center" one without k offset
    dynamical_matrix!(H, sswt, zero(Vec3); branch=2)
    δE₁ = -real(tr(view(H, 1:L, 1:L))) / 2Natoms
    # @show δE₁

    for branch in 1:3
        dynamical_matrix!(H, sswt, zero(Vec3); branch)
        @show -real(tr(view(H, 1:L, 1:L))) / 2Natoms
    end

    # Integrate zero-point energy over 𝐪 ∈ [0, 1]³ in reshaped RLU, and sum
    # over 3 spiral branches (q-k, q, q+k)
    δE₂ = hcubature((0,0,0), (1,1,1); opts...) do q_reshaped
        total = 0.0
        for branch in 1:3
            dynamical_matrix!(H, sswt, q_reshaped; branch)
            ωs = bogoliubov!(V, H)
            total += sum(view(ωs, 1:L))
        end
        return total / (3 * 2Natoms)  # normalize over atoms and branches
    end

    @show δE₂[1]

    # Error bars in δE₂[2] are discarded
    return δE₁ + δE₂[1]
end

function magnetization_lswt_correction_dipole(sswt::SpinWaveTheorySpiral; opts...)
    L = nbands(sswt.swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    δS = hcubature((0,0,0), (1,1,1); opts...) do q
        dynamical_matrix!(H, sswt, Vec3(q); branch=2)
        bogoliubov!(V, H)
        return SVector{L}(-norm2(view(V, L+i, 1:L)) for i in 1:L)
    end

    # Error bars in δS[2] are discarded
    return δS[1]
end
