module PlottingExt

using Sunny
import Sunny: Mat3, Vec3, orig_crystal, natoms
using LinearAlgebra
import Statistics
import Makie
import Brillouin

include("PlottingUtils.jl")
include("ViewCrystal.jl")
include("PlotSpins.jl")
include("PlotIntensities.jl")
include("ViewQSpace.jl")

function __init__()
    if Sunny.is_pkg_loaded(:WGLMakie)
        @info """
        WGLMakie (the web backend) is experimental. If you encounter graphics
        problems, try loading GLMakie instead of WGLMakie from a fresh session.
        Issue tracker: https://github.com/SunnySuite/Sunny.jl/issues/211.
        """
    end
end

end
