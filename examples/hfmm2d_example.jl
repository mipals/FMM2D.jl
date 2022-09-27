using FMM2D
using LinearAlgebra

# Simple example for the FMM2D Library
eps = 10.0^(-5)  # Tolerance
zk  = 1.1 + im*0 # Wavenumber

# Source-to-source,
n = 200
sources = rand(2,n)
charges = rand(ComplexF64,n)
vals = hfmm2d(eps,zk,sources,charges=charges,pg=1)
vals.pot

# Source-to-target
ntargets = 300
targets  = rand(2,ntargets)
vals = hfmm2d(eps,zk,sources;charges=charges,targets=targets,pgt=1)
vals.pottarg
