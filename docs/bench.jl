# To benchmark example builds: `julia --project=@. bench.jl`

import IOCapture, Sunny, GLMakie

examples = filter(endswith(".jl"), readdir(pkgdir(Sunny, "examples"), join=true))

for desc in ("Cold start", "Warm start")
    println(desc)
    for ex in examples
        print(rpad(basename(ex), 20), ": ")
        @time IOCapture.capture() do
            include(ex)
        end
    end
    println()
end
