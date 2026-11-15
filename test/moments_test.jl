@testset "Exact spectral traces" begin
    s = small_system()
    (; H, Ja, Jb, Jc, psi, T, NH, NC) = s
    exact = real.([tr(t) / NH for t in T])
    @test KPMru.kpm_1d(H, NC, NH; psi_in=psi) ≈ exact atol=1e-12
    per_vector = KPMru.kpm_1d(H, NC, NH; psi_in=psi, avg_output=false)
    @test per_vector ≈ hcat(diag.(T)...) atol=1e-12
    @test KPMru.kpm_1d_current(H, Ja, NC, NH; psi_in=psi) ≈
          real.([tr(Ja * t) / NH for t in T]) atol=1e-12
    current_samples = KPMru.kpm_1d_current(H, Ja, NC, NH; psi_in=psi, avg_output=false)
    @test current_samples ≈ hcat([diag(Ja * t) for t in T]...) atol=1e-12
    for site in 1:NH
        @test KPMru.ldos_mu(H, NC, site) ≈ real.([t[site, site] for t in T]) atol=1e-12
    end

    # Index order follows the operators acting on the right and left vectors.
    exact2 = [tr(T[m] * Jb * T[n] * Ja) / NH for n in 1:NC, m in 1:NC]
    for arr_size in (3, 4)
        mu2 = KPMru.kpm_2d(H, Ja, Jb, NC, NH, NH; psi_in=psi, arr_size=arr_size)
        @test mu2 ≈ exact2 atol=1e-12
    end
    for (parity, predicate) in ((:EVEN, iseven), (:ODD, isodd))
        mu2 = KPMru.kpm_2d(H, Ja, Jb, NC, NH, NH; psi_in=psi, moment_parity=parity)
        @test mu2 ≈ [predicate(n-m) ? exact2[n,m] : 0 for n in 1:NC, m in 1:NC] atol=1e-12
    end
    @test_throws ArgumentError KPMru.kpm_2d(H, Ja, Jb, NC, NH, NH;
                                          psi_in=psi, moment_parity=:INVALID)
    @test_throws AssertionError KPMru.kpm_1d(H, NC, NH; psi_in=psi[:, 1:2])
    @test_throws AssertionError KPMru.kpm_1d(H, NC, NH; psi_in=psi, psi_in_l=psi)
    @test_throws AssertionError KPMru.kpm_1d(H, NC, NH; psi_in_l=psi)

    NC3 = 4
    mu3 = KPMru.kpm_3d(H, Ja, Jb, Jc, NC3, NH, NH; psi_in=psi)
    exact3 = [tr(T[m] * Jc * T[n] * Jb * T[p] * Ja) / NH
              for n in 1:NC3, m in 1:NC3, p in 1:NC3]
    @test mu3 ≈ exact3 atol=1e-12
end
