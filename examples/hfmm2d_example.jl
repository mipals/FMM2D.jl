using FMM2D

# Simple example for the FMM2D Library
thresh = 10.0^(-5)          # Tolerance
zk     = rand(ComplexF64)   # Wavenumber

# Source-to-source,
n = 200
sources = rand(2,n)
charges = rand(ComplexF64,n)
vals    = hfmm2d(thresh,zk,sources,charges=charges,pg=1)

# Source-to-target
ntargets = 300
targets  = rand(2,ntargets)
vals     = hfmm2d(thresh,zk,sources;charges=charges,targets=targets,pgt=1)
vals.pottarg
