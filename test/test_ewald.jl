println("test_ewald")


using Parameters

# Many of these tests will be broken and need re-tooling to new interface

function test_compression(A, Acomp)
    ndim = div(ndims(A), 2) - 1
    nb = size(A, 1)
    latsize = size(A)[2:1+ndim]
    for i in CartesianIndices(latsize)
        for j in CartesianIndices(latsize)
            for ib in 1:nb
                for jb in 1:nb
                    if !(A[ib, i, jb, j] ≈ Acomp[ib, jb, i - j])
                        return false
                    end
                end
            end
        end
    end

    return true
end

function test_ewald_NaCl()
    lat_vecs = [1.0 0   0;
                0   1.0 0;
                0   0   1.0]
    b_vecs = [zeros(3)]
    latsize = [2, 2, 2]
    lattice = Sunny.Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1, -1, 1, -1, 1, 1, -1], (1, 2, 2, 2))

    ewald_result = Sunny.ewald_sum_monopole(lattice, sys.sites, extent=30)
    direct_result = Sunny.direct_sum_monopole(lattice, sys.sites, extent=30)

    answer = -1.7475645946331822
    @test isapprox(answer, direct_result / 4; rtol=1e-7)
    @test isapprox(answer, ewald_result / 4; rtol=1e-7)
end

function test_ewald_CsCl()
    lat_vecs = [1.0 0   0;
                0   1.0 0;
                0   0   1.0]
    b_vecs = [[0.,  0.,  0.],
              [0.5, 0.5, 0.5]]
    latsize = [1, 1, 1]
    lattice = Sunny.Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1], (2, 1, 1, 1))

    ewald_result = Sunny.ewald_sum_monopole(lattice, sys.sites, extent=30)
    # direct_result = Sunny.direct_sum_monopole(lattice, sys.sites, extent=30)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/4)
    # direct_result *= √(3/4)

    answer = -1.76267477307099
    # @test isapprox(answer, direct_result; rtol=1e-7)
    @test isapprox(answer, ewald_result; rtol=1e-7)
end

function test_ewald_ZnS()
    lat_vecs = [0.0 0.5 0.5;
                  0.5 0.0 0.5;
                  0.5 0.5 0.0]
    b_vecs = [[0.,  0.,  0.],
              [0.25, 0.25, 0.25]]
    latsize = [1, 1, 1]
    lattice = Sunny.Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, -1], (2, 1, 1, 1))

    ewald_result = Sunny.ewald_sum_monopole(sys, extent=30)
    # direct_result = direct_sum_monopole(sys, extent=200)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= √(3/16)
    # direct_result *= √(3/16)

    answer = -1.63805505338879
    # @test isapprox(answer, direct_result; rtol=1e-7)
    @test isapprox(answer, ewald_result; rtol=1e-7)
end

function test_ewald_ZnSB4()
    a = 1.
    c = √(8/3) * a
    u = 3/8

    lat_vecs = [ 0.5a    0.5a    0.0;
                -0.5*√3a 0.5*√3a 0.0;
                 0.0     0.0       c]
    b_vecs = [[0.5a,  0.5/√3*a,  0.],
              [0.5a, -0.5/√3*a, 0.5c,],
              [0.5a, 0.5/√3*a, u*c],
              [0.5a, -0.5/√3*a, (0.5+u)*c]]
    b_vecs = [[1/3, 2/3, 0],
              [2/3, 1/3, 0.5],
              [1/3, 2/3, 3/8],
              [2/3, 1/3, 7/8]]
    latsize = [1, 1, 1]
    lattice = Sunny.Lattice(lat_vecs, b_vecs, latsize)
    sys = ChargeSystem(lattice)
    sys.sites .= reshape([1, 1, -1, -1], (4, 1, 1, 1))

    ewald_result = Sunny.ewald_sum_monopole(lattice, sys.sites, extent=30)
    # direct_result = direct_sum_monopole(sys, extent=200)

    # Madelung constants are reported relative to the
    #  nearest-neighbor distance in the crystal
    ewald_result *= u * c
    # direct_result *= u * c

    answer = -1.64132162737
    @test isapprox(answer, ewald_result / 2; rtol=1e-7)
end

"Self-energy of a physical dipole with moment p, and displacement d=2ϵ"
function _dipole_self_energy(; p::Float64=1.0, ϵ::Float64=0.1)
    d, q = 2ϵ, p/2ϵ
    return -q^2 / d
end

"""Approximates a dipolar `SpinSystem` by generating a monopolar `ChargeSystem` consisting of
    opposite charges Q = ±1/(2ϵ) separated by displacements d = 2ϵp centered on the original
    lattice sites.

    I think this function may be wrong now.
"""
function _approx_dip_as_mono(sys::SpinSystem{D, L, Db}; ϵ::Float64=0.1) :: ChargeSystem{D, L, Db} where {D, L, Db}
    @unpack sites, lattice = sys

    # Need to expand the underlying unit cell to the entire system size
    new_lat_vecs = lattice.size' .* lattice.lat_vecs
    new_latsize = @SVector ones(Int, D)

    frac_transform = inv(new_lat_vecs)

    new_nbasis = 2 * prod(size(sites))
    new_sites = zeros(new_nbasis, 1, 1, 1)
    new_basis = Vector{SVector{D, Float64}}()
    sizehint!(new_basis, new_nbasis)

    ib = 1
    for idx in eachindex(lattice)
        @inbounds r = lattice[idx]
        @inbounds p = sites[idx]

        # Add new charges as additional basis vectors
        push!(new_basis, frac_transform * SVector{3}(r .+ ϵ * p))
        push!(new_basis, frac_transform * SVector{3}(r .- ϵ * p))

        # Set these charges to ±1/2ϵ
        new_sites[ib, 1, 1, 1]   =  1 / (2ϵ)
        new_sites[ib+1, 1, 1, 1] = -1 / (2ϵ)

        ib += 2
    end

    new_basis = map(b->mod.(b, 1), new_basis)

    new_lattice = Sunny.Lattice(new_lat_vecs, new_basis, new_latsize)

    return ChargeSystem{D, L, Db}(new_lattice, new_sites)
end

"""
Tests that `ewald_sum_monopole` and `ewald_sum_dipole` give consistent results
 up to approximation error.
TODO: This version uses existing functions, but suffers from catastrophic
 cancellation issues due to the self-energy exploding. Should write a 
 separate monopole ewald which skips these self-energy interactions.
"""
function test_mono_dip_consistent()
    lat_vecs = [1.0 0   0;
                0   1.0 0;
                0   0   1.0]
    b_vecs = [[0.,  0.,  0.],
              [0.5, 0.5, 0.5]]
    latsize = [1, 1, 1]
    cryst = Crystal(lat_vecs, b_vecs)
    sys = SpinSystem(cryst, Sunny.Interaction[], latsize)
    rand!(sys)

    dip_ewald = Sunny.ewald_sum_dipole(sys.lattice, sys.sites; extent=50)

    csys = _approx_dip_as_mono(sys; ϵ=0.001)
    mono_ewald = Sunny.ewald_sum_monopole(csys.lattice, csys.sites; extent=50)
    dip_self_en = _dipole_self_energy(; ϵ=0.001)

    println("Dipole Ewald Energy: $(dip_ewald)")
    println("Monopole Ewald Energy: $(mono_ewald - dip_self_en)")
    @assert isapprox(dip_ewald, mono_ewald - dip_self_en; rtol=1e-4)
end

# Beck claims this should be independent of η, and converge to -2.837297
# I find that it only does for large η?
function test_ξ_sum(;extent=2, η=1.0) :: Float64
    extent_ixs = CartesianIndices(ntuple(_->-extent:extent, Val(3)))

    real_space_sum = 0.0
    recip_space_sum = 0.0

    for JKL in extent_ixs
        JKL = convert(SVector, JKL)
        k = 2π * JKL

        if all(JKL .== 0)
            continue
        end

        dist = norm(JKL)
        kdist = 2π * dist

        real_space_sum += erfc(η * dist) / dist
        recip_space_sum += exp(-kdist^2 / (4η^2)) / kdist^2
    end

    return real_space_sum + 4π * recip_space_sum - 2η/√π - π/η^2
end