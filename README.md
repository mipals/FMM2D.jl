# FMM2D.jl

[![Build Status](https://github.com/mipals/FMM2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/FMM2D.jl/actions/workflows/CI.yml?query=branch%3Amain)

FMM2D.jl is a Julia interface for computing N-body interactions using the [Flatiron Institute's FMM2D library](https://github.com/flatironinstitute/fmm2d/).

Currently, the interface only supports Helmholtz problems.

## Helmholtz FMM
Let $c_j \in \mathbb{C},\ j = 1,\dots,N$ denote a collection of charge strengths, $v_j \in \mathbb{C},\ j = 1,\dots,N$ denote a collection of dipole strengths, and $\mathbf{d}_j\in\mathbb{R}^2,\ j = 1,\dots,N$ denote the corresponding dipole orientation vectors. Furthermore, $k \in \mathbb{C}$ denote the wave number.

The Fast Multipole Method is $O(n)$ algorithm that can be used to approximate the Helmholtz potential $u$ (and its gradient and Hessian) caused by the presence of a collection of sources $\mathbf{x}_j$ at target position $\mathbf{x}$. The direct computation of the potential is 

$$
u(\mathbf{x}) = \sum_{j=1}^{N} c_jH_0^{(1)}(k\|\mathbf{x} - \mathbf{x}_j\|) - v_j\mathbf{d}_j\cdot\nabla H_0^{(1)}(k\|\mathbf{x} - \mathbf{x}_j\|),
$$

where $H_0^{(1)}$ is the Hankel function of the first kind of order 0. When $\mathbf{x} = \mathbf{x}_j$ the $j$ th term is dropped from the sum.


## Example
```julia
using FMM2D

# Simple example for the FMM2D Library
thresh = 10.0^(-5)          # Tolerance
zk     = rand(ComplexF64)   # Wavenumber

# Source-to-source,
n = 200
sources = rand(2,n)
charges = rand(ComplexF64,n)
vals = hfmm2d(thresh,zk,sources,charges=charges,pg=1)
```
