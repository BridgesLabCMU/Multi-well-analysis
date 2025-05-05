#!/usr/bin/env julia
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "Pipeline.jl"))
using .Pipeline

function julia_main()
	try
        Pipeline.pipeline()
        return 0 # if things finished successfully
	catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
	end
end
exit(julia_main())
