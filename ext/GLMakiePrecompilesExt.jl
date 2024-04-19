module GLMakiePrecompilesExt

import Sunny, GLMakie

import PrecompileTools as PT
PT.@setup_workload begin
    PT.@compile_workload begin
        cryst = Sunny.bcc_crystal()
        display(Sunny.view_crystal(cryst); visible=false)

        fig, ax, plt = GLMakie.meshscatter([GLMakie.Point3f(0,0,0), GLMakie.Point3f(1,0,0)])
        inspector = GLMakie.DataInspector(ax; font="Deja Vu Sans Mono")
        GLMakie.Makie.show_data(inspector, plt, 1)
        display(fig; visible=false)
        
        GLMakie.closeall()
    end
end

end
