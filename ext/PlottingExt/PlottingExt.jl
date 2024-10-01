module PlottingExt

using Sunny
import Sunny: Mat3, Vec3, orig_crystal, natoms
using LinearAlgebra
import Statistics
import Makie

let warned = false
    global warn_wglmakie() = begin
        if !warned && string(Makie.current_backend()) == "WGLMakie"
            @info """
            Using the WGLMakie graphics backend. If you encounter graphics problems,
            try restarting the Julia session and load GLMakie instead of WGLMakie.
            Issue tracker: https://github.com/SunnySuite/Sunny.jl/issues/211.
            """
        end
        warned = true
    end
end

include("PlottingUtils.jl")
include("ViewCrystal.jl")
include("PlotSpins.jl")
include("PlotIntensities.jl")
include("ViewQSpace.jl")

end
