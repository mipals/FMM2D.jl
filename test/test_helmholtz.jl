using FMM2D
using SpecialFunctions
using LinearAlgebra
using Random
using Test

Random.seed!(0)

@testset "Helmholtz FMM 2D (complex)" begin
    zk = rand(ComplexF64)
    N = 1000
    sources = rand(2, N)
    charges = rand(N) + 1im*rand(N)
    dipstrs = rand(N) + 1im*rand(N)
    dipvecs = ones(2, N)/sqrt(2)
    # Direct computations (N^2)
    refpot = FMM2D.h2ddir(zk=zk,sources=sources,targets=sources,dipvecs=dipvecs,dipstr=dipstrs,charges=charges,pgt=3)

    M = 100
    targets = rand(2, M)
    # Direct computations (N^2)
    refpottarg = FMM2D.h2ddir(zk=zk,sources=sources,targets=targets,dipvecs=dipvecs,dipstr=dipstrs,charges=charges,pgt=3)
    for tol_exp=-14:-1
        tolerance = 0.49*10.0^tol_exp
        U = hfmm2d(eps=tolerance,zk=zk,sources=sources, charges=charges, dipstr=dipstrs, dipvecs=dipvecs,pg=3)
        pot_relerr  = norm(U.pot  - refpot.pottarg)  / norm(refpot.pottarg)
        grad_relerr = norm(U.grad - refpot.gradtarg) / norm(refpot.gradtarg)
        hess_relerr = norm(U.hess - refpot.hesstarg) / norm(refpot.hesstarg)
        @test pot_relerr  < tolerance
        @test grad_relerr < tolerance
        @test hess_relerr < tolerance

        Utarg = hfmm2d(eps=tolerance,zk=zk,sources=sources, charges=charges, dipstr=dipstrs, dipvecs=dipvecs,pgt=3,targets=targets)
        pot_relerrtarg  = norm(Utarg.pottarg  - refpottarg.pottarg)  / norm(refpottarg.pottarg)
        grad_relerrtarg = norm(Utarg.gradtarg - refpottarg.gradtarg) / norm(refpottarg.gradtarg)
        hess_relerrtarg = norm(Utarg.hesstarg - refpottarg.hesstarg) / norm(refpottarg.hesstarg)
        @test pot_relerrtarg  < tolerance
        @test grad_relerrtarg < tolerance
        @test hess_relerrtarg < tolerance
    end
end
