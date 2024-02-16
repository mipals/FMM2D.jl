using SafeTestsets

@safetestset "Aqua testing      " begin include("test_aqua.jl")      end
@safetestset "Laplace Wrappers  " begin include("test_laplace.jl")   end
@safetestset "Helmholtz Wrappers" begin include("test_helmholtz.jl") end
# @safetestset "Stokes Wrappers   " begin include("test_stokes.jl")    end
