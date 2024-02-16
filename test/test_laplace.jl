using FMM2D
using LinearAlgebra
using Random
using Test

Random.seed!(0)

function direct_self(source, charges, dipvecs, dipstrs)
    if size(charges,2) == 1
        charges = transpose(charges)
    end
    if size(dipstrs,2) == 1
        dipstrs = transpose(dipstrs)
    end
    nd,n = size(charges)
    pot = zeros(eltype(charges),nd,n)
    N = size(source, 2)
    for k=1:nd
        charge = charges[k,:]
        dipvec = dipvecs[2k-1:2k,:]
        dipstr = dipstrs[k,:]
        for i=1:N
            for j=1:N
                if i == j
                    continue
                end
                r1 = source[1,j]-source[1,i]
                r2 = source[2,j]-source[2,i]
                rsq = r1^2 + r2^2
                rdotq = r1*dipvec[1,j] + r2*dipvec[2,j]
                pot[k,i] += charge[j]*log(rsq)/2 + dipstr[j]*rdotq/rsq
            end
        end
    end
    if nd == 1
        return transpose(pot)
    else
        return pot
    end
end
function direct_targ(source, charges, dipvecs, dipstrs, target)
    if size(charges,2) == 1
        charges = transpose(charges)
    end
    if size(dipstrs,2) == 1
        dipstrs = transpose(dipstrs)
    end
    nd = size(charges,1)
    M = size(target, 2)
    N = size(source, 2)
    pot = zeros(eltype(charges), nd, M)
    for k in 1:nd
        charge = charges[k,:]
        dipvec = dipvecs[2k-1:2k,:]
        dipstr = dipstrs[k,:]
        for i=1:M
            for j=1:N
                r1 = source[1,j] - target[1,i]
                r2 = source[2,j] - target[2,i]
                rsq = r1^2 + r2^2
                rdotq = r1*dipvec[1,j] + r2*dipvec[2,j]
                pot[k,i] += charge[j]*log(rsq)/2 + dipstr[j]*rdotq/rsq
            end
        end
    end
    if nd == 1
        return transpose(pot)
    else
        return pot
    end
end

@testset "Laplace FMM 2D (complex)" begin
    N = 1000
    sources = rand(2, N)
    charges = rand(ComplexF64,N)
    dipstrs = rand(ComplexF64,N)
    dipvecs = ones(2, N)/sqrt(2)
    # Direct computations (N^2)
    refpot = direct_self(sources,charges,dipvecs,dipstrs)

    M = 100
    targets = rand(2, M)
    # Direct computations (N^2)
    refpottarg = direct_targ(sources,charges,dipvecs,dipstrs,targets)
    for tol_exp=-14:-1
        tolerance = 0.49*10.0^tol_exp
        U = lfmm2d(eps=tolerance,sources=sources,charges=charges, dipstr=dipstrs, dipvecs=dipvecs,pg=2)
        pot_relerr  = norm(U.pot  - refpot)  / norm(refpot)
        # grad_relerr = norm(U.grad - refpot.gradtarg) / norm(refpot.gradtarg)
        # hess_relerr = norm(U.hess - refpot.hesstarg) / norm(refpot.hesstarg)
        @test pot_relerr  < tolerance
        # @test grad_relerr < tolerance
        # @test hess_relerr < tolerance

        Utarg = lfmm2d(eps=tolerance,sources=sources, charges=charges, dipstr=dipstrs, dipvecs=dipvecs,targets=targets,pgt=1)
        pot_relerrtarg  = norm(Utarg.pottarg  - refpottarg)  / norm(refpottarg)
        # grad_relerrtarg = norm(Utarg.gradtarg - refpottarg.gradtarg) / norm(refpottarg.gradtarg)
        # hess_relerrtarg = norm(Utarg.hesstarg - refpottarg.hesstarg) / norm(refpottarg.hesstarg)
        @test pot_relerrtarg  < tolerance
        # @test grad_relerrtarg < tolerance
        # @test hess_relerrtarg < tolerance
    end
end

@testset "Laplace FMM 2D (real)" begin
    N = 1000
    sources = rand(2, N)
    charges = rand(N)
    dipstrs = rand(N)
    dipvecs = ones(2, N)/sqrt(2)
    # Direct computations (N^2)
    refpot = direct_self(sources,charges,dipvecs,dipstrs)

    M = 100
    targets = rand(2, M)
    # Direct computations (N^2)
    refpottarg = direct_targ(sources,charges,dipvecs,dipstrs,targets)
    for tol_exp=-14:-1
        tolerance = 0.49*10.0^tol_exp
        U = rfmm2d(eps=tolerance,sources=sources,charges=charges, dipstr=dipstrs, dipvecs=dipvecs,pg=2)
        pot_relerr  = norm(U.pot  - refpot)  / norm(refpot)
        # grad_relerr = norm(U.grad - refpot.gradtarg) / norm(refpot.gradtarg)
        # hess_relerr = norm(U.hess - refpot.hesstarg) / norm(refpot.hesstarg)
        @test pot_relerr  < tolerance
        # @test grad_relerr < tolerance
        # @test hess_relerr < tolerance

        Utarg = rfmm2d(eps=tolerance,sources=sources, charges=charges, dipstr=dipstrs, dipvecs=dipvecs,targets=targets,pgt=1)
        pot_relerrtarg  = norm(Utarg.pottarg  - refpottarg)  / norm(refpottarg)
        # grad_relerrtarg = norm(Utarg.gradtarg - refpottarg.gradtarg) / norm(refpottarg.gradtarg)
        # hess_relerrtarg = norm(Utarg.hesstarg - refpottarg.hesstarg) / norm(refpottarg.hesstarg)
        @test pot_relerrtarg  < tolerance
        # @test grad_relerrtarg < tolerance
        # @test hess_relerrtarg < tolerance
    end
end
