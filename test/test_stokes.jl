using FMM2D
using LinearAlgebra
using Random
using Test

Random.seed!(0)

# Brute-force O(N^2) reference, straight from the kernel definitions documented in
# `stfmm2d` (the library returns the raw sum, WITHOUT the 1/2pi scaling):
#   G_ij(x,y) = r_i r_j / (2 r^2) - delta_ij log(r) / 2
#   T_ijk(x,y) = -2 r_i r_j r_k / r^4
#   u_i(x) = sum_m G_ij sigma_j + sum_m T_ijk mu_j nu_k
function stokes_direct(sources, stoklet, strslet, strsvec, targets)
    n = size(sources, 2)
    nt = size(targets, 2)
    pot = zeros(2, nt)
    for t in 1:nt
        x = @view targets[:, t]
        for m in 1:n
            y = @view sources[:, m]
            r = x - y
            rr = dot(r, r)
            rr == 0 && continue # drop self-interaction
            d = sqrt(rr)
            for i in 1:2, j in 1:2
                G = r[i] * r[j] / (2rr) - (i == j) * log(d) / 2
                pot[i, t] += G * stoklet[j, m]
                for k in 1:2
                    T = -2 * r[i] * r[j] * r[k] / rr^2
                    pot[i, t] += T * strslet[j, m] * strsvec[k, m]
                end
            end
        end
    end
    return pot
end

@testset "Stokes FMM 2D" begin
    N = 1000
    sources = rand(2, N)
    stoklet = rand(2, N)
    strslet = rand(2, N)
    strsvec = rand(2, N)

    M = 100
    targets = rand(2, M) .+ 2.0 # well separated from the sources

    # reference for source->target and source->source evaluations
    ref_targ = stokes_direct(sources, stoklet, strslet, strsvec, targets)
    ref_self = stokes_direct(sources, stoklet, strslet, strsvec, sources)

    for tol_exp in -12:-1
        tol = 0.49 * 10.0^tol_exp

        # Stokeslet + stresslet, source -> target
        U = stfmm2d(;
            eps = tol, sources = sources,
            stoklet = stoklet, strslet = strslet, strsvec = strsvec,
            targets = targets, ppregt = 1,
        )
        @test norm(U.pottarg - ref_targ) / norm(ref_targ) < tol

        # Stokeslet + stresslet, source -> source (self term dropped)
        Us = stfmm2d(;
            eps = tol, sources = sources,
            stoklet = stoklet, strslet = strslet, strsvec = strsvec,
            ppreg = 1,
        )
        @test norm(Us.pot - ref_self) / norm(ref_self) < tol
    end

    # Stokeslet only and stresslet only, to exercise each piece independently
    tol = 1.0e-12
    ref_stok = stokes_direct(sources, stoklet, zeros(2, N), zeros(2, N), targets)
    Ustok = stfmm2d(; eps = tol, sources = sources, stoklet = stoklet, targets = targets, ppregt = 1)
    @test norm(Ustok.pottarg - ref_stok) / norm(ref_stok) < tol

    ref_strs = stokes_direct(sources, zeros(2, N), strslet, strsvec, targets)
    Ustrs = stfmm2d(;
        eps = tol, sources = sources,
        strslet = strslet, strsvec = strsvec, targets = targets, ppregt = 1,
    )
    @test norm(Ustrs.pottarg - ref_strs) / norm(ref_strs) < tol
end
