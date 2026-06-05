"""
```julia
    anyfail, n, nt, ifstoklet, ifstrslet = (
        stfmm2dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd) )
```

Input checking for the Stokes routine. Checks the sizes of the input arrays and
the compatibility of the various flags, `nd`, and provided arrays. If something
is off, a warning is issued and `anyfail` is set to `true`.

Output:
* `anyfail`   - boolean, `true` if any checks fail, `false` otherwise
* `n`         - integer, number of sources
* `nt`        - integer, number of targets
* `ifstoklet` - integer, `1` if there are Stokeslet strengths, `0` otherwise
* `ifstrslet` - integer, `1` if there are stresslet strengths, `0` otherwise
"""
function stfmm2dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd)
    anyfail = false

    if (size(sources,1) != 2)
        @warn "sources array is wrong shape, computing nothing"
        anyfail = true
    end
    if (nd < 0)
        @warn "nd has invalid value, computing nothing"
        anyfail = true
    end
    if (ppreg > 3 || ppreg < 0)
        @warn "flag ppreg not in expected range, computing nothing"
        anyfail = true
    end
    if (ppregt > 3 || ppregt < 0)
        @warn "flag ppregt not in expected range, computing nothing"
        anyfail = true
    end

    n = div(length(sources),2)
    nt = 0

    if targets !== nothing
        if (size(targets,1) != 2)
            @warn "targets array is wrong shape, computing nothing"
            anyfail = true
        end
        nt = div(length(targets),2)
    end

    if stoklet !== nothing
        if (div(length(stoklet),nd) != 2*n)
            @warn "size of stoklet array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifstoklet = 1
    else
        ifstoklet = 0
    end

    if (strslet !== nothing || strsvec !== nothing)
        if (strsvec === nothing || strslet === nothing)
            @warn "if stresslets in calculation, need both strslet and strsvec arrays, computing nothing"
            anyfail = true
        end
        if (div(length(strslet),nd) != 2*n)
            @warn "size of strslet array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        if (div(length(strsvec),nd) != 2*n)
            @warn "size of strsvec array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifstrslet = 1
    else
        ifstrslet = 0
    end

    if ifstoklet == 0 && ifstrslet == 0
        @warn "no Stokeslets or stresslets provided, doing nothing"
        anyfail = true
    end

    if (ppreg != 1 && ppreg != 2 && ppreg != 3) && (ppregt != 1 && ppregt != 2 && ppregt != 3)
        @warn "no output requested, doing nothing"
        anyfail = true
    end

    if (ppregt == 1 || ppregt == 2 || ppregt == 3) && targets === nothing
        @warn "target values requested but no targets provided, proceeding anyway"
    end

    return anyfail, n, nt, ifstoklet, ifstrslet
end

"""
```julia
    vals = stfmm2d(eps,sources;
            stoklet=nothing,strslet=nothing,strsvec=nothing,
            targets=nothing,ppreg=0,ppregt=0,nd=1)
```

This subroutine computes the N-body Stokes interactions in two dimensions. For a
source ``y`` and target ``x``, let ``r_i = x_i - y_i`` and ``r = \\|x-y\\|``. The
Stokeslet ``G_{ij}`` and its associated pressure tensor ``P_j`` (without the
``1/2\\pi`` scaling) are
```math
    G_{ij}(x,y) = \\frac{r_i r_j}{2 r^2} - \\frac{\\delta_{ij}}{2}\\log(r), \\qquad
    P_j(x,y) = \\frac{r_j}{r^2},
```
and the (Type I) stresslet ``T_{ijk}`` and its associated pressure tensor
``\\Pi_{jk}`` (without the ``1/2\\pi`` scaling) are
```math
    T_{ijk}(x,y) = -2\\frac{r_i r_j r_k}{r^4}, \\qquad
    \\Pi_{jk}(x,y) = -\\frac{\\delta_{jk}}{r^2} + 2\\frac{r_j r_k}{r^4}.
```
The velocity ``u`` and pressure ``p`` are then
```math
    u_i(x) = \\sum_{m} G_{ij}(x,y^{(m)})\\,\\sigma^{(m)}_j
           + \\sum_{m} T_{ijk}(x,y^{(m)})\\,\\mu^{(m)}_j\\,\\nu^{(m)}_k,
```
```math
    p(x) = \\sum_{m} P_j(x,y^{(m)})\\,\\sigma^{(m)}_j
         + \\sum_{m} \\Pi_{jk}(x,y^{(m)})\\,\\mu^{(m)}_j\\,\\nu^{(m)}_k,
```
where ``\\sigma^{(m)}`` are the Stokeslet strengths (`stoklet`), ``\\mu^{(m)}``
the stresslet strengths (`strslet`) and ``\\nu^{(m)}`` the stresslet orientation
vectors (`strsvec`), located at the sources ``y^{(m)}``. When ``x = y^{(m)}`` the
``m`` th term is dropped from the sum.

# Input
* `eps::Float64` precision requested
* `sources::Array{Float64}` size (2,n) source locations (``y^{(m)}``)

# Optional Keyword Input
* `stoklet::Array{Float64}` size (nd,2,n) or (2,n) Stokeslet strengths (``\\sigma``)
* `strslet::Array{Float64}` size (nd,2,n) or (2,n) stresslet strengths (``\\mu``)
* `strsvec::Array{Float64}` size (nd,2,n) or (2,n) stresslet orientations (``\\nu``)
* `targets::Array{Float64}` size (2,nt) target locations (``x``)
* `ppreg::Integer` source eval flag.
    + Do not compute any quantities at sources if `ppreg == 0`
    + Velocity (``u``) at sources evaluated if `ppreg == 1`.
    + Velocity and pressure (``p``) at sources evaluated if `ppreg == 2`
    + Velocity, pressure, and gradient (``\\nabla u``) at sources evaluated if `ppreg == 3`
* `ppregt::Integer` target eval flag (same meaning as `ppreg` but at targets).
* `nd::Integer` number of densities

# Output
`vals::FMMVals` with the fields (sizes given for `nd == 1`)
* `vals.pot::Array{Float64}` size (2,n) velocity at source locations if requested
* `vals.pre::Array{Float64}` size (n) pressure at source locations if requested
* `vals.grad::Array{Float64}` size (2,2,n) velocity gradient at sources if requested
* `vals.pottarg::Array{Float64}` size (2,nt) velocity at target locations if requested
* `vals.pretarg::Array{Float64}` size (nt) pressure at target locations if requested
* `vals.gradtarg::Array{Float64}` size (2,2,nt) velocity gradient at targets if requested
* `vals.ier` error flag as returned by FMM2D library; `0` indicates a successful call.
"""
function stfmm2d(;eps::Float64,sources::Array{Float64},
        stoklet::TFN=nothing,strslet::TFN=nothing,strsvec::TFN=nothing,
        targets::TFN=nothing,ppreg::Integer=0,ppregt::Integer=0,nd::Integer=1)

    # default values
    vals = FMMVals()

    zero  = Float64(0)
    zero2 = zeros(Float64,2)

    pot      = zero
    pre      = zero
    grad     = zero
    pottarg  = zero
    pretarg  = zero
    gradtarg = zero

    # check inputs
    anyfail, n, nt, ifstoklet, ifstrslet = (
        stfmm2dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd))

    if anyfail
        return vals
    end

    if (ifstoklet == 0); stoklet = zero2 end
    if (ifstrslet == 0); strslet = zero2; strsvec = zero2 end
    if (nt == 0);        targets = zero2 end

    # allocate memory for return values
    if ppreg == 1 || ppreg == 2 || ppreg == 3
        pot = (nd > 1) ? zeros(Float64,nd,2,n) : zeros(Float64,2,n)
    end
    if ppreg == 2 || ppreg == 3
        pre = (nd > 1) ? zeros(Float64,nd,n) : zeros(Float64,n)
    end
    if ppreg == 3
        grad = (nd > 1) ? zeros(Float64,nd,2,2,n) : zeros(Float64,2,2,n)
    end
    if ppregt == 1 || ppregt == 2 || ppregt == 3
        pottarg = (nd > 1) ? zeros(Float64,nd,2,nt) : zeros(Float64,2,nt)
    end
    if ppregt == 2 || ppregt == 3
        pretarg = (nd > 1) ? zeros(Float64,nd,nt) : zeros(Float64,nt)
    end
    if ppregt == 3
        gradtarg = (nd > 1) ? zeros(Float64,nd,2,2,nt) : zeros(Float64,2,2,nt)
    end

    ier = Integer(0)

    # Calling the Fortran function
    # subroutine stfmm2d(nd, eps, nsource, source,
    #     ifstoklet, stoklet, ifstrslet, strslet, strsvec,
    #     ifppreg, pot, pre, grad, ntarg, targ,
    #     ifppregtarg, pottarg, pretarg, gradtarg, ier)
    ccall((:stfmm2d_,libfmm2d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fd,Fi,
                                      Fd,Fd,Fi,Fd,Fd,Fd,Fi,
                                      Fd,Fi,Fd,Fd,Fd,Fi),
          nd,eps,n,sources,ifstoklet,stoklet,ifstrslet,
          strslet,strsvec,ppreg,pot,pre,grad,nt,targets,ppregt,
          pottarg,pretarg,gradtarg,ier)

    if (ier != 0); @warn "libfmm2d had an error, see vals.ier"; end
    vals.ier = ier

    # Save computed values to output struct
    if ppreg == 1 || ppreg == 2 || ppreg == 3; vals.pot = pot end
    if ppreg == 2 || ppreg == 3; vals.pre = pre end
    if ppreg == 3; vals.grad = grad end
    if ppregt == 1 || ppregt == 2 || ppregt == 3; vals.pottarg = pottarg end
    if ppregt == 2 || ppregt == 3; vals.pretarg = pretarg end
    if ppregt == 3; vals.gradtarg = gradtarg end

    return vals
end
