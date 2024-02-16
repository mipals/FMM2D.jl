"""
```julia
    vals = lfmm2d(eps,zk,sources;
            charges=nothing,dipstr=nothing,dipvecs=nothing,targets=nothing,pg=0,pgt=0,nd=1)
```
This subroutine computes the N-body Laplace interactions in two dimensions where the
interaction kernel is given by log(r) and its gradients.
```math
    u(x) = \\sum_{j=1}^{N} c_{j} * log\\(\\|x-x_{j}\\|\\) + d_{j}v_{j} \\cdot \\nabla( log(\\|x-x_{j}\\|) ),
```

where   ``c_{j}`` are the charge densities,
        ``d_{j}`` are the dipole strength vectors
        ``v_{j}`` are the dipole orientation vectors,
        ``x_{j}`` are the source locations.
when ``x=x_{j}``, the jth term is dropped from the sum

# Input
* `eps::Float64` precision requested
* `sources::Array{Float64}` size (2,n) source locations (``x_{j}``)

# Optional Keyword Input
* `charges::Array{ComplexF64}` size (nd,n)
* `dipstr::Array{ComplexF64}`  size (nd,n)
* `dipvecs::Array{Float64}` size (nd,2,n) or (2,n) dipole orientation vectors (``d_{j}``)
* `targets::Array{Float64}` size (2,nt) target locations (``x``)
* `pg::Integer` source eval flag.
    + Do not compute any quantities at sources if `pg == 0`
    + Potential (``u``) at sources evaluated if `pg == 1`.
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
    + Potential, gradient, and Hessian (``\\nabla \\nabla u``) at sources evaluated if `pg == 3`
* `pgt::Integer` target eval flag.
    + Do not compute any quantities at targets if `pgt == 0`
    + Potential at targets evaluated if `pgt == 1`.
    + Potential and gradient at targets evaluated if `pgt == 2`
    + Potential, gradient, and Hessian at targets evaluated if `pgt == 3`
* `nd::Integer` number of densities
Note: if all default values are used for optional input, nothing is computed.

# Output
`vals::FMMVals` with the fields
Note that the Hessian is returned as a length 4 vector at each point with the second derivatives
in the order: ``\\partial_{xx}``, ``\\partial_{yy}``, ``\\partial_{xy}``, ``\\partial_{yx}``.
* `vals.pot::Array{ComplexF64}` size (nd,n) or (n) potential at source locations if requested
* `vals.grad::Array{ComplexF64}` size (nd,2,n) or (2,n) gradient at source locations if requested
* `vals.hess::Array{ComplexF64}` size (nd,4,n) or (4,n) Hessian at source locations if requested
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,2,nt) or (2,nt) gradient at target locations if requested
* `vals.hesstarg::Array{ComplexF64}` size (nd,4,nt) or (4,nt) Hessian at target locations if requested
* `vals.ier` error flag as returned by FMM2D library. A value of 0 indicates a successful call.
Non-zero values may indicate insufficient memory available. See the documentation for the FMM2D library.
If not set (`nothing`), then FMM2D library was never called.
"""
function lfmm2d(;eps::Float64,sources::Array{Float64},
        charges::TCN=nothing,dipvecs::TFN=nothing,dipstr::TCN=nothing,targets::TFN=nothing,
        pg::Integer=0,pgt::Integer=0,nd::Integer=1)


    # default values
    vals = FMMVals()

    zero = ComplexF64(0)

    #
    pot = zero
    grad = zero
    hess = zero
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm2dinputcheck(sources,charges,dipvecs,dipstr,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if dipstr === nothing && dipvecs !== nothing
        dipstr = ones(eltype(charges),n)
    end

    if (ifcharge == 0); charges = zero  end
    if (ifdipole == 0); dipvecs = 0.0   end
    if (ifdipole == 0); dipstr  = zero  end
    if (nt == 0);       targets = 0.0   end

    # allocate memory for return values

    if pg == 1 || pg == 2 || pg == 3
        if nd > 1
            pot = zeros(ComplexF64,nd,n)
        else
            pot = zeros(ComplexF64,n)
        end
    end

    if pg == 2 || pg == 3
        if nd > 1
            grad = zeros(ComplexF64,nd,2,n)
        else
            grad = zeros(ComplexF64,2,n)
        end
    end

    if pg == 3
        if nd > 1
            hess = zeros(ComplexF64,nd,3,n)
        else
            hess = zeros(ComplexF64,3,n)
        end
    end

    if pgt == 1 || pgt == 2 || pgt == 3
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2 || pgt == 3
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,2,nt)
        else
            gradtarg = zeros(ComplexF64,2,nt)
        end
    end

    if pgt == 3
        if nd > 1
            hesstarg = zeros(ComplexF64,nd,3,nt)
        else
            hesstarg = zeros(ComplexF64,3,nt)
        end
    end


    ier = Integer(0)
    iper = Integer(0)

    # Calling the Fortran function
    # subroutine lfmm2d(nd,eps,ns,sources,ifcharge,charge,
    #                   ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
    #                   nt,targ,ifpghtarg,pottarg,gradtarg,
    #                   hesstarg,ier)
    ccall((:lfmm2d_,libfmm2d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fc,Fi,
                                        Fc,Fd,Fi,Fi,Fc,Fc,Fc,Fi,
                                        Fd,Fi,Fc,Fc,Fc,Fi),
        nd,eps,n,sources,ifcharge,charges,ifdipole,
        dipstr,dipvecs,iper,pg,pot,grad,hess,nt,targets,pgt,
        pottarg,gradtarg,hesstarg,ier)

    if (vals.ier != 0); @warn "libfmm2d had an error, see vals.ier"; end

    # Save computed values to output struct
    if pg == 1 || pg == 2 || pg == 3; vals.pot = pot end
    if pg == 2 || pg == 3; vals.grad = grad end
    if pg == 3; vals.hess = hess end
    if pgt == 1 || pgt == 2 || pgt == 3; vals.pottarg = pottarg end
    if pgt == 2 || pgt == 3; vals.gradtarg = gradtarg end
    if pgt == 3; vals.hesstarg = hesstarg end

    return vals

end


"""
```julia
    vals = rfmm2d(eps,zk,sources;
            charges=nothing,dipstr=nothing,dipvecs=nothing,targets=nothing,pg=0,pgt=0,nd=1)
```
This subroutine computes the N-body Laplace interactions in two dimensions where the
interaction kernel is given by log(r) and its gradients.
```math
    u(x) = \\sum_{j=1}^{N} c_{j} * log\\(\\|x-x_{j}\\|\\) + d_{j}v_{j} \\cdot \\nabla( log(\\|x-x_{j}\\|) ),
```

where   ``c_{j}`` are the charge densities,
        ``d_{j}`` are the dipole strength vectors
        ``v_{j}`` are the dipole orientation vectors,
        ``x_{j}`` are the source locations.
when ``x=x_{j}``, the jth term is dropped from the sum

# Input
* `eps::Float64` precision requested
* `sources::Array{Float64}` size (2,n) source locations (``x_{j}``)

# Optional Keyword Input
* `charges::Array{ComplexF64}` size (nd,n)
* `dipstr::Array{ComplexF64}`  size (nd,n)
* `dipvecs::Array{Float64}` size (nd,2,n) or (2,n) dipole orientation vectors (``d_{j}``)
* `targets::Array{Float64}` size (2,nt) target locations (``x``)
* `pg::Integer` source eval flag.
    + Do not compute any quantities at sources if `pg == 0`
    + Potential (``u``) at sources evaluated if `pg == 1`.
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
    + Potential, gradient, and Hessian (``\\nabla \\nabla u``) at sources evaluated if `pg == 3`
* `pgt::Integer` target eval flag.
    + Do not compute any quantities at targets if `pgt == 0`
    + Potential at targets evaluated if `pgt == 1`.
    + Potential and gradient at targets evaluated if `pgt == 2`
    + Potential, gradient, and Hessian at targets evaluated if `pgt == 3`
* `nd::Integer` number of densities
Note: if all default values are used for optional input, nothing is computed.

# Output
`vals::FMMVals` with the fields
Note that the Hessian is returned as a length 4 vector at each point with the second derivatives
in the order: ``\\partial_{xx}``, ``\\partial_{yy}``, ``\\partial_{xy}``, ``\\partial_{yx}``.
* `vals.pot::Array{ComplexF64}` size (nd,n) or (n) potential at source locations if requested
* `vals.grad::Array{ComplexF64}` size (nd,2,n) or (2,n) gradient at source locations if requested
* `vals.hess::Array{ComplexF64}` size (nd,4,n) or (4,n) Hessian at source locations if requested
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,2,nt) or (2,nt) gradient at target locations if requested
* `vals.hesstarg::Array{ComplexF64}` size (nd,4,nt) or (4,nt) Hessian at target locations if requested
* `vals.ier` error flag as returned by FMM2D library. A value of 0 indicates a successful call.
Non-zero values may indicate insufficient memory available. See the documentation for the FMM2D library.
If not set (`nothing`), then FMM2D library was never called.
"""
function rfmm2d(;eps::Float64,sources::Array{Float64},
        charges::TFN=nothing,dipvecs::TFN=nothing,dipstr::TFN=nothing,targets::TFN=nothing,
        pg::Integer=0,pgt::Integer=0,nd::Integer=1)


    # default values
    vals = FMMVals()

    zero = Float64(0)

    #
    pot = zero
    grad = zero
    hess = zero
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm2dinputcheck(sources,charges,dipvecs,dipstr,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if dipstr === nothing && dipvecs !== nothing
        dipstr = ones(n)
    end

    if (ifcharge == 0); charges = zero  end
    if (ifdipole == 0); dipvecs = 0.0   end
    if (ifdipole == 0); dipstr  = zero  end
    if (nt == 0);       targets = 0.0   end

    # allocate memory for return values

    if pg == 1 || pg == 2 || pg == 3
        if nd > 1
            pot = zeros(nd,n)
        else
            pot = zeros(n)
        end
    end

    if pg == 2 || pg == 3
        if nd > 1
            grad = zeros(nd,2,n)
        else
            grad = zeros(2,n)
        end
    end

    if pg == 3
        if nd > 1
            hess = zeros(nd,4,n)
        else
            hess = zeros(4,n)
        end
    end

    if pgt == 1 || pgt == 2 || pgt == 3
        if nd > 1
            pottarg = zeros(nd,nt)
        else
            pottarg = zeros(nt)
        end
    end

    if pgt == 2 || pgt == 3
        if nd > 1
            gradtarg = zeros(nd,2,nt)
        else
            gradtarg = zeros(2,nt)
        end
    end

    if pgt == 3
        if nd > 1
            hesstarg = zeros(nd,4,nt)
        else
            hesstarg = zeros(4,nt)
        end
    end


    ier = Integer(0)
    iper = Integer(0)

    # Calling the Fortran function
    # subroutine rfmm2d(nd,eps,ns,sources,ifcharge,charge,
    #               ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
    #               nt,targ,ifpghtarg,pottarg,gradtarg,
    #               hesstarg,ier)
    ccall((:rfmm2d_,libfmm2d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fd,Fi,
                                    Fd,Fd,Fi,Fi,Fd,Fd,Fd,Fi,
                                    Fd,Fi,Fd,Fd,Fd,Fi),
          nd,eps,n,sources,ifcharge,charges,ifdipole,
          dipstr,dipvecs,iper,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg,ier)

    if (vals.ier != 0); @warn "libfmm2d had an error, see vals.ier"; end

    # Save computed values to output struct
    if pg == 1 || pg == 2 || pg == 3; vals.pot = pot end
    if pg == 2 || pg == 3; vals.grad = grad end
    if pg == 3; vals.hess = hess end
    if pgt == 1 || pgt == 2 || pgt == 3; vals.pottarg = pottarg end
    if pgt == 2 || pgt == 3; vals.gradtarg = gradtarg end
    if pgt == 3; vals.hesstarg = hesstarg end

    return vals

end
