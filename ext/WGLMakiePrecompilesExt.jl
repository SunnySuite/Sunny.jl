module WGLMakiePrecompilesExt

# Disabling, because it does not seem to meaningfully reduce latency.

#=
import Sunny, WGLMakie

import PrecompileTools as PT
PT.@setup_workload begin
    PT.@compile_workload begin
        cryst = Sunny.bcc_crystal()
        try
            WGLMakie.@compile Sunny.view_crystal(cryst)
        catch e
            println(e)
        end
    end
end
=#

end
