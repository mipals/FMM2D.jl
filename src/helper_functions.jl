
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
