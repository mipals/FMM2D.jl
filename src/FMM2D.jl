module FMM2D

using FMM2D_jll

export hfmm2d

# fortran input/return types
Fd = Ref{Float64}
Fi = Ref{Int32}
Fc = Ref{ComplexF64}

# common input types
TFN = Union{Array{Float64},Nothing}
TCN = Union{Array{ComplexF64},Nothing}

mutable struct FMMVals
    pot
    grad
    hess
    pottarg
    gradtarg
    hesstarg
    ier
end

function FMMVals()
    return FMMVals(nothing,nothing,nothing,nothing,nothing,nothing,0)
end


"""
```julia
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm2dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd) )
```

Input checking for Helmholtz and Laplace routines.
Checks the sizes of these arrays and the compatibility
of the various flags, nd, and provided arrays.
If something is off, a warning is issued and
`anyfail` is set to true. Non-fatal mistakes result in
a warning but `anyfail` remains false.

Output:
* `anyfail` - boolean is true if any checks fail, false otherwise
* n - integer, number of sources
* nt - integer, number of targets
* ifcharge - integer, is 1 if there are charges, 0 otherwise
* ifdipole - integer, is 1 if there are dipoles, 0 otherwise
"""
function scalarfmm2dinputcheck(sources,charges,dipvecs,dipstr,targets,pg,pgt,nd)
    anyfail = false

    if (size(sources,1) != 2)
        @warn "sources array is wrong shape, computing nothing"
        anyfail = true
    end
    if (nd < 0)
        @warn "nd has invalid value, computing nothing"
        anyfail = true
    end
    if (pg > 3 || pg < 0)
        @warn "flag pg not in expected range, computing nothing"
        anyfail = true
    end
    if (pgt > 3 || pgt < 0)
        @warn "flag pgt not in expected range, computing nothing"
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

    if charges !== nothing
        if (div(length(charges),nd) != n)
            @warn "size of charges array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifcharge = 1
    else
        ifcharge = 0
    end

    if dipvecs !== nothing
        if (div(length(dipvecs),nd) != n*2)
            @warn "size of dipvecs array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        if dipstr !== nothing
            if (div(length(dipstr),nd) != n)
                @warn "size of dipstr array incompatible with sources array and nd parameter, computing nothing"
                anyfail = true
            end
        end
        ifdipole = 1
    else
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        anyfail = true
    end

    if (pg != 1 && pg != 2 && pg != 3) && (pgt != 1 && pgt != 2 && pgt != 3)
        @warn "no output requested, doing nothing"
        anyfail = true
    end

    if (pgt == 1 || pgt == 2 || pgt == 3) && targets === nothing
        @warn "target values requested but no targets provided, proceeding anyway"
    end

    return anyfail, n, nt, ifcharge, ifdipole
end



function hfmm2d(eps::Float64,zk::Union{Float64,ComplexF64},
                sources::Array{Float64};
                charges::TCN=nothing,dipvecs::TFN=nothing,dipstr::TCN=nothing,
                targets::TFN=nothing,pg::Integer=0,pgt::Integer=0,
                nd::Integer=1)

    # casting zk to a complex number
    zk = complex(zk)

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
        dipstr = ones(ComplexF64,n)
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
            hess = zeros(ComplexF64,nd,4,n)
        else
            hess = zeros(ComplexF64,4,n)
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
            hesstarg = zeros(ComplexF64,nd,4,nt)
        else
            hesstarg = zeros(ComplexF64,4,nt)
        end
    end

    # actually call the function
    # hfmm2d(nd,eps,zk,ns,sources,ifcharge,
    #          charge, ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,nt,
    #           targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
    # dipstr = dipole strength?

    ier = Integer(0)
    iper = Integer(0)

    ccall((:hfmm2d_,libfmm2d),Cvoid,(Fi,Fd,Fc,Fi,Fd,Fi,Fc,Fi,
                                    Fc,Fd,Fi,Fi,Fc,Fc,Fc,Fi,
                                    Fd,Fi,Fc,Fc,Fc,Fi),
          nd,eps,zk,n,sources,ifcharge,charges,ifdipole,
          dipstr,dipvecs,iper,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg,ier)

    if (vals.ier != 0); @warn "libfmm2d had an error, see vals.ier"; end

    # load requested values
    if pg == 1 || pg == 2 || pg == 3; vals.pot = pot end
    if pg == 2 || pg == 3; vals.grad = grad end
    if pg == 3; vals.hess = hess end
    if pgt == 1 || pgt == 2 || pgt == 3; vals.pottarg = pottarg end
    if pgt == 2 || pgt == 3; vals.gradtarg = gradtarg end
    if pgt == 3; vals.hesstarg = hesstarg end

    return vals

end


function h2ddir(zk::Union{ComplexF64,Float64},sources::Array{Float64},
                targets::Array{Float64};
                charges::TCN=nothing,dipvecs::TFN=nothing,dipstr::TCN=nothing,
                pgt::Integer=0,nd::Integer=1,
                thresh::Float64=1e-15)

    zk = complex(zk)

    # default values
    vals = FMMVals()

    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)

    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    pg = 0
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm2dinputcheck(sources,charges,dipvecs,dipstr,targets,pg,pgt,nd))

    if dipstr === nothing && dipvecs !== nothing
        dipstr = ones(ComplexF64,n)
    end

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero  end
    if (ifdipole == 0); dipvecs = 0.0   end
    if (ifdipole == 0); dipstr  = zero  end
    if (nt == 0); targets = 0.0 end

    # allocate memory for return values

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
            hesstarg = zeros(ComplexF64,nd,4,nt)
        else
            hesstarg = zeros(ComplexF64,4,nt)
        end
    end

    # hfmm2dpart_direct(nd,istart,iend,jstart,jend,zk,source,ifcharge,
    #       charge,ifdipole,dipstr,dipvec,targ,ifpgh,
    #       pot,grad,hess,thresh)
    #     This subroutine adds the contribuition due to sources
    #     istart to iend in the source array at the expansion centers
    #     jstart to jend in the target array to the computed velocities
    #     and gradients. Note that contributions for sources
    #     within thresh of the targets are not added to the potential
    # dispatch to appropriate wrapper

    ccall((:hfmm2dpart_direct_,libfmm2d),Cvoid,(Fi,Fi,Fi,Fi,Fi,Fc,Fd,Fi,
                                            Fc,Fi,Fc,Fd,Fd,Fi,
                                            Fc,Fc,Fc,Fd),
                            nd,1,n,1,nt,zk,sources,ifcharge,
                            charges,ifdipole,dipstr,dipvecs,targets,pgt,
                            pottarg,gradtarg,hesstarg,thresh)


    if pgt == 1 || pgt == 2 || pgt == 3; vals.pottarg = pottarg end
    if pgt == 2 || pgt == 3; vals.gradtarg = gradtarg end
    if pgt == 3; vals.hesstarg = hesstarg end

    return vals

end

end
