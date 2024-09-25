
@testitem "Parse nxs" begin
    @test Sunny.parse_long_name("[0.5H,0.3K,0.1H]") == [0.5, 0.3, 0.1]
    # TODO: Add stripped down .nxs files of various flavors
end


@testitem "Binning" begin
    using LinearAlgebra

    # Setup an arbitrary qgrid
    latvecs = lattice_vectors(1, 1, 10, 90, 90, 90)
    cryst = Crystal(latvecs, [[0, 0, 0]])
    qgrid = q_space_grid(cryst, [0, 1.2, 0.5], range(0, 1, 100), [0, -8, 0.1], (-3,-1))

    dqx = diff(qgrid.qs, dims=1)[1]
    dqy = diff(qgrid.qs, dims=2)[1]

    @test_throws "in-plane" Sunny.specify_transverse_binning(qgrid, dqx, 0.1)

    params = Sunny.specify_transverse_binning(qgrid, dqx × dqy, 0.1)

    # There should be exactly one transverse bin
    @test params.numbins[3] == 1

    # Energy axis should be fully integrated by default
    @test isinf(params.binwidth[4])

    # Test that moving by one gridpoint in the QGrid moves by one binwidth in
    # binning coordinate space
    @test (params.covectors[1:3,1:3] * dqx)[1] ≈ params.binwidth[1]
    @test (params.covectors[1:3,1:3] * dqy)[2] ≈ params.binwidth[2]

    # Make a system with only a few momenta available
    sys = System(cryst, [1 => Moment(s=1/2, g=2)], :dipole; dims=(5,5,1))
    randomize_spins!(sys)

    static_corrs = SampledCorrelationsStatic(sys; measure=ssf_trace(sys))
    add_sample!(static_corrs,sys)

    res_interp = intensities_static(static_corrs, qgrid)
    res_bin = Sunny.binned_intensities(static_corrs, params)

    # The binning and interpolating results should be very different because
    # only a few momenta are available. The binned data will be very sparse,
    # while interpolation makes it dense
    @test count(iszero.(res_bin.data)) > 10 * count(iszero.(res_interp.data))

    # Returning to a better behaved q-grid, we should get agreement
    unit_params = Sunny.unit_resolution_binning_parameters(static_corrs)
    res_bin_unit = Sunny.binned_intensities(static_corrs, unit_params)

    # Construct the equivalent grid of q points 
    bcs = Sunny.axes_bincenters(unit_params)
    unit_qgrid_qpts = Sunny.QPoints([[qx,qy,qz] for qx = bcs[1], qy = bcs[2], qz = bcs[3]][:])
    res_interp_unit = intensities_static(static_corrs,unit_qgrid_qpts)

    ratio = res_bin_unit.data[:,:,1,1] ./ reshape(res_interp_unit.data,5,5)

    # The results should be proportional
    @test sum((ratio .- ratio[1]).^2) < 1e-12

    # The constant of proportionality should be the bin volume.
    @test ratio[1] ≈ prod(unit_params.binwidth[1:3])

    # This is the expected sum rule for this setup (TODO: expression in terms of
    # N, S, etc)
    @test sum(res_bin_unit.data) ≈ 1
end
