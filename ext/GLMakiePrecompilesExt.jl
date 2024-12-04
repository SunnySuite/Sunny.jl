module GLMakiePrecompilesExt

import Sunny, Makie, GLMakie
import PrecompileTools as PT

# Julia 1.11.1 broke extension precompiles, but this was fixed in 1.11.2.
# https://github.com/JuliaLang/julia/issues/56204#issuecomment-2439588043
@static if VERSION != v"1.11.1"

PT.@setup_workload begin
    PT.@compile_workload begin
        cryst = Sunny.bcc_crystal()
        GLMakie.display(Sunny.view_crystal(cryst); visible=false)

        fig, ax, plt = GLMakie.meshscatter([GLMakie.Point3f(0,0,0), GLMakie.Point3f(1,0,0)])
        inspector = GLMakie.DataInspector(ax; font="Deja Vu Sans Mono")
        GLMakie.Makie.show_data(inspector, plt, 1)
        GLMakie.display(fig; visible=false)

        GLMakie.closeall()
    end
end

end

end
