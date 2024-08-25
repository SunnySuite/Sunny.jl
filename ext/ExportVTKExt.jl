module ExportVTKExt

using Sunny
using LinearAlgebra
using StaticArrays
import WriteVTK 


function Sunny.export_vtk(filename,sys; coordinates = :physical, log_scale = false)
    if !ispath(filename)
        mkdir(filename)
    end
    # In SU(N) mode, `sys.dipoles` are expectation values
    dipole_data_name = sys.mode == :SUN ? "expectation_value_dipole" : "dipole"

    Ni, Nj, Nk = sys.dims
    latvecs = if coordinates == :physical
        sys.crystal.latvecs
    elseif coordinates == :lattice
        1. * I(3)
    end
    x = [dot([i,j,k], latvecs[1,:]) for i = 1:Ni, j = 1:Nj, k = 1:Nk]
    y = [dot([i,j,k], latvecs[2,:]) for i = 1:Ni, j = 1:Nj, k = 1:Nk]
    z = [dot([i,j,k], latvecs[3,:]) for i = 1:Ni, j = 1:Nj, k = 1:Nk]

    # SQTODO: Dipoles point in an un-physical direction in :lattice coordinate mode
    dipole_data = [[sys.dipoles[i,j,k,l] for i = 1:Ni, j = 1:Nj, k = 1:Nk] for l = 1:natoms(sys.crystal)]

    ## Save dipole vectors

    saved_files = WriteVTK.vtk_multiblock(joinpath(filename,"dipoles")) do vtm
        for l in 1:length(sys.crystal.positions)
            offset = latvecs * sys.crystal.positions[l]
            WriteVTK.vtk_grid(vtm,x .+ offset[1], y .+ offset[2], z .+ offset[3]) do vtk
                vtk[dipole_data_name] = dipole_data[l]
                vtk["latvecs",WriteVTK.VTKFieldData()] = latvecs
            end
        end
    end

    ## Save S(Q) correlations
    
    ic = instant_correlations(sys)
    add_sample!(ic, sys)
    params = unit_resolution_binning_parameters(ic)
    formula = intensity_formula(ic, :trace)

    # Counts is always just one because of unit_resolution
    signal, _ = intensities_binned(ic,params; formula)

    if log_scale
        signal .= log10.(signal)
    end

    saved_files_more = export_vtk(joinpath(filename,"correlations"), params, signal; dims_kept = [1,2,3])
    append!(saved_files, saved_files_more)
end

"""
    export_vtk(filename,params::BinningParameters,data)

Export a VTK-compatible file to `filename` (do not include file extension when
specifying the file name) which contains the `data` as VTK Cell Data on a grid
parameterized by `params`.

At least one axis of the [`BinningParameters`](@ref) must be integrated over,
since VTK does not support 4D data. See `integrate_axes!`.
"""
function Sunny.export_vtk(filename,params::BinningParameters,data;dims_kept = nothing)
    # Storing the bin *edges* as the grid in the VTK file
    # so that the data is treated as cell data (correct for histogram data)
    edges = Vector{AbstractRange{Float64}}(undef,0)

    dims_integrated_over, dims_kept = begin
      if isnothing(dims_kept) # Auto-detect which dimensions are singleton
          dims_integrated_over = Int64[]
          for dim in 1:4
              if params.numbins[dim] > 1
                  push!(edges, Sunny.axes_binedges(params.binstart[dim],params.binend[dim],params.binwidth[dim])[1])
              else
                  push!(dims_integrated_over, dim)
              end
          end
          dims_integrated_over, setdiff(1:4,dims_integrated_over)
      else
          for dim in dims_kept
              push!(edges, Sunny.axes_binedges(params.binstart[dim],params.binend[dim],params.binwidth[dim])[1])
          end
          setdiff(1:4,dims_kept), dims_kept
      end
    end

    if length(dims_integrated_over) == 0
        error("VTK doesn't support 4D data")
    end

    # Remove the integrated dimensions from the data
    data = dropdims(data;dims = Tuple(dims_integrated_over))

    nd = length(dims_kept)
    bcs = axes_bincenters(params)[dims_kept]
    centers = [SVector{nd,Float64}([bcs[i][ix[i]] for i = 1:nd]) for ix in CartesianIndices(size(data))]

    WriteVTK.vtk_grid(filename, edges...) do vtk
        vtk["data"] = data
        vtk["bin_center"] = centers
    end
end

end
