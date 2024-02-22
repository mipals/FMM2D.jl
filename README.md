# FMM2D.jl

[![Build Status](https://github.com/mipals/FMM2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mipals/FMM2D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mipals/FMM2D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mipals/FMM2D.jl)

FMM2D.jl is a Julia interface for computing N-body interactions using the [Flatiron Institute's FMM2D library](https://github.com/flatironinstitute/fmm2d/).

Currently, the wrapper only wraps the Helmholtz and Laplace functionalities.

## Helmholtz FMM
Let $c_s \in \mathbb{C},\ s = 1,\dots,N$ denote a collection of charge strengths, $v_s \in \mathbb{C},\ s = 1,\dots,N$ denote a collection of dipole strengths, and $d_s\in\mathbb{R}^2,\ s = 1,\dots,N$ denote the corresponding dipole orientation vectors. Furthermore, $k \in \mathbb{C}$ denotes the wave number. The Helmholtz potential $u$ caused by the presence of a collection of $M$ sources ($x_s$) at $N$ target positions ($x_t$) is computed as

$$
    u\left(x_t\right) = \sum_{s=1}^{M} c_sH_0^{(1)}(k\|x_t - x_s\|) - v_sd_s\cdot\nabla H_0^{(1)}(k\|x_t - x_s\|), \quad t = 1,\dots, N
$$

where $H_0^{(1)}$ is the Hankel function of the first kind of order 0. When $x = x_j$ the $j$ th term is dropped from the sum. Performing this summation would scale as $O(NM)$, but using the Flatiron Insitutes Fast Multipole Library a linear scaling of $O((N + M)\text{log}(\varepsilon^{-1}))$ can be achieved with $\varepsilon$ being the desired relative precision. Note that the library also includes the option for computing the gradient and Hessian of the potential.


### Example
```julia
using FMM2D

# Simple example for the FMM2D Library
thresh = 10.0^(-5)          # Tolerance
zk     = rand(ComplexF64)   # Wavenumber

# Source-to-source
n = 200
sources = rand(2,n)
charges = rand(ComplexF64,n)
pg = 3 # Evaluate potential, gradient, and Hessian at the sources
vals = hfmm2d(eps=thresh,zk=zk,sources=sources,charges=charges,pg=pg)
vals.pot
vals.grad
vals.hess

# Source-to-target
m = 200
targets = rand(2,m)
pgt = 3 # Evaluate potential, gradient, and Hessian at the targets
vals = hfmm2d(targets=targets,eps=thresh,zk=zk,sources=sources,charges=charges,pgt=pgt)
vals.pottarg
vals.gradtarg
vals.hesstarg
```

## Laplace
The Laplace problem in 2D have the following form

$$
    u(x) = \sum_{j=1}^{N} \left[c_{j} \text{log}\left(\|x-x_{j}\|\right) - d_{j}v_{j} \cdot \nabla( \text{log}(\|x-x_{j}\|) )\right],
$$

In the case of complex charges and dipole strengths ($c_j, v_j \in \mathbb{C}^n$) the function call `lfmm2d` has to be used. In the case of real charges and dipole strengths ($c_j, v_j \in \mathbb{R}^n$) the function call `rfmm2d` has to be used.


### Example
```julia
using FMM2D

# Simple example for the FMM2D Library
thresh = 10.0^(-5)          # Tolerance

# Source-to-source
n = 200
sources = rand(2,n)
charges = rand(ComplexF64,n)
dipvecs = randn(2,n)
dipstr = rand(ComplexF64,n)
pg = 3 # Evaluate potential, gradient, and Hessian at the sources
vals = lfmm2d(eps=thresh,sources=sources,charges=charges,dipvecs=dipvecs,dipstr=dipstr,pg=pg)
vals.pot
vals.grad
vals.hess

# Source-to-target
m = 100
targets = rand(2,m)
pgt = 3 # Evaluate potential, gradient, and Hessian at the targets
vals = lfmm2d(targets=targets,eps=thresh,sources=sources,charges=charges,dipvecs=dipvecs,dipstr=dipstr,pgt=pgt)
vals.pottarg
vals.gradtarg
vals.hesstarg
```


<!-- In addition, the `FMM2D` library also includes the following sum

$$
    u(z) = \sum_{j=1}^{N} \left[c_{j} \text{log}\left(z-\varepsilon_j\right) - \frac{v_j}{z - \varepsilon_j}\right],
$$

where $c_j \ in$ -->

<!-- ## Biharmonic

Let $c_j = (c_{j,1}, c_{j,2}) \in \mathbb{C}^2,\ j=1,2,\dots, N$ denote a collection of charge strengths and $v_j = (v_{j,1}, v_{j,2}, v_{j,3}) \in \mathbb{C}^3, j=1,2,\dots,N$ denote a collection of dipole strengths. 

The Biharmonic FMM computes the potential $u(x)$ 

$$
    u(z) = \sum_{j=1}^N \left[c_{j,1}\text{log}\left(\|z - \varepsilon_j\|\right) + c_{j,2}\frac{z - \varepsilon_j}{\overline{z - \varepsilon_j}} + \frac{v_{j,1}}{z - \varepsilon_j} + \frac{v_{j,3}}{\overline{z - \varepsilon_j}} + v_{j,2}\frac{z - \varepsilon_j}{\left(\overline{z - \varepsilon_j}\right)^2}\right],
$$

as well as its gradient $(P_z\frac{\mathrm{d}}{\mathrm{d}z}, P_{\overline{z}}\frac{\mathrm{d}}{\mathrm{d}z}, \frac{\mathrm{d}}{\mathrm{d}\overline{z}})$ given by

$$
\begin{aligned}
    P_z\frac{\mathrm{d}}{\mathrm{d}z}u(z) &= \sum_{j=1}^{N}\left[\frac{c_{j,1}}{z - \varepsilon_j} - \frac{v_{j,1}}{\left(z - \varepsilon_j\right)^2}\right]\\
    P_{\overline{z}}\frac{\mathrm{d}}{\mathrm{d}z}u(z) &= \sum_{j=1}^{N}\left[\frac{c_{j,2}}{\overline{z - \varepsilon_j}} - \frac{v_{j,2}}{\left(\overline{z - \varepsilon_j}\right)^2}\right]\\
    \frac{\mathrm{d}}{\mathrm{d}\overline{z}}u(z) &= \sum_{j=1}^{N}\left[\frac{c_{j,1}}{\overline{z - \varepsilon_j}} - c_{j,2}\frac{z - \varepsilon_{j}}{\left(\overline{z - \varepsilon_j}\right)^2}- v_{j,3}\frac{z - \varepsilon_{j}}{\left(\overline{z - \varepsilon_j}\right)^2} - 2v_{j,2}\frac{z - \varepsilon_{j}}{\left(\overline{z - \varepsilon_j}\right)^3}\right]
\end{aligned}.
$$ 

at the source and target locations. When $z = \varepsilon_j$ the $j$ th term is dropped from the sum. Note here that $z = x_1 + i x_2$ and $\varepsilon_j = x_{j,1} + ix_{j,2}$. -->



## Stokes

$$
    u(x) = \sum_{j=1}^NG^\text{stok}(x,x_j)c_j + d_j\cdot T^\text{stok}(x,x_j)\cdot v_j
$$

$$
    p(x) = \sum_{j=1}^NP^\text{stok}(x,x_j)c_j + d_j\cdot \Pi^\text{stok}(x,x_j)\cdot v_j^\top
$$

## Related Package
[FMMLIB2D.jl](https://github.com/ludvigak/FMMLIB2D.jl) interfaces the [FMMLIB2D](https://github.com/zgimbutas/fmmlib2d) library which the [FMM2D library improves on](https://fmm2d.readthedocs.io/en/latest/). 

