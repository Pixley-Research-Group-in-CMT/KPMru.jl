@testset "DOS reconstruction and derivatives" begin
    # mu_0=1, other moments zero is the exactly known arcsine distribution.
    mu = [1.0; zeros(15)]
    a = 2.0
    grid = [-2.2, -0.7, 0.0, 0.4, 2.2]
    E, rho = KPMru.dos(mu, a; E_grid=grid, N_tilde=length(grid))
    expected = [abs(x) < a ? 1 / (pi * sqrt(a^2-x^2)) : 0.0 for x in grid]
    @test E == grid
    @test rho ≈ expected
    _, drho = KPMru.dos(mu, a; E_grid=grid, N_tilde=length(grid), dE_order=1)
    @test drho ≈ [abs(x) < a ? x / (pi * (a^2-x^2)^(3/2)) : 0.0 for x in grid]
    @test KPMru.dos0(mu, a) ≈ 1/(pi*a)
    @test KPMru.dos0(mu, a; dE_order=2) ≈ 1/(pi*a^3)
    @test KPMru.dos0(mu, a; NC=8) ≈ KPMru.dos0(mu[1:8], a)
    @test_throws AssertionError KPMru.dos(mu, -a)
    @test_throws AssertionError KPMru.dos(mu, a; E_grid=grid, E_range=[-1.0,1.0])
    @test_throws AssertionError KPMru.dos(mu, a; E_grid=grid, N_tilde=2)
    @test_throws AssertionError KPMru.dos(mu, a; dE_order=2)

    # A nontrivial finite Hamiltonian: exact moments, Jackson positivity and mass.
    s = small_system(NC=32)
    moments = KPMru.kpm_1d(s.H, s.NC, s.NH; psi_in=s.psi)
    theta = ((1:512) .- 0.5) .* (pi/512)
    energies = cos.(theta)
    _, density = KPMru.dos(moments, 1.0; E_grid=energies, N_tilde=length(energies))
    @test minimum(density) >= -1e-12
    # Integrate with E=cos(theta), avoiding the endpoint square-root singularity.
    @test sum(density .* sin.(theta)) * pi/512 ≈ 1 atol=1e-12
    _, at_zero = KPMru.dos(moments, 1.0; E_grid=[0.0], N_tilde=1)
    @test only(at_zero) ≈ KPMru.dos0(moments, 1.0)
    h = 1e-4
    _, near_zero = KPMru.dos(moments, 1.0; E_grid=[-h, 0.0, h], N_tilde=3)
    @test (near_zero[1]-2near_zero[2]+near_zero[3])/h^2 ≈
          KPMru.dos0(moments, 1.0; dE_order=2) rtol=1e-5
end
