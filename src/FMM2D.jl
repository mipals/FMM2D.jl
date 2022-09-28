module FMM2D

# Importing binaries
using FMM2D_jll

# Exporting interface
export hfmm2d

# Fortran input/return types
Fd = Ref{Float64}
Fi = Ref{Int32}
Fc = Ref{ComplexF64}

# Common input types
TFN = Union{Array{Float64},Nothing}
TCN = Union{Array{ComplexF64},Nothing}

# Return struct
mutable struct FMMVals
    pot         # Potential at sources
    grad        # Gradient at souces
    hess        # Hessian at sources
    pottarg     # Potential at targets
    gradtarg    # Gradient at targets
    hesstarg    # Hessian at targets
    ier         # Error-indicator
end
function FMMVals()
    return FMMVals(nothing,nothing,nothing,nothing,nothing,nothing,0)
end

# Include helper functions
include("helper_functions.jl")
# Include wrappers
include("helmholtz_wrappers.jl")

end
