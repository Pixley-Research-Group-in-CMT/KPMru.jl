using JLD2
using Test
using KPMru
using LinearAlgebra
using Random

@load joinpath(@__DIR__, "test_data", "ops_small.jld2") H_norm a Jx
Random.seed!(402)

mu = KPMru.kpm_1d(H_norm, 2048, 21)
rho_0 = KPMru.dos0(mu, a)
d2rho_0 = KPMru.dos0(mu, a; dE_order=2)
println("rho_0 = $(rho_0)")
println("d2rho_0 = $(d2rho_0)")
@test abs(rho_0) < 1e-4


NC=4
NR=2
NH=H_norm.n

Gamma = zeros(ComplexF64, NC, NC, NC)

psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPMru.normalize_by_col(psi_in_l, NR; centering=false)
psi_in_r = psi_in_l

KPMru.kpm_3d!(H_norm, Jx, Jx, Jx, NC, NR, NH, Gamma, psi_in_l, psi_in_r)
ω = 0.1
# Exceptions from the actual integration call must fail this test.
cpge_val = KPMru.cpge(Gamma, NC, ω)
@test isfinite(cpge_val)
@test isfinite(d2rho_0)
@test KPMru.cpge(2Gamma, NC, ω) ≈ 2cpge_val
