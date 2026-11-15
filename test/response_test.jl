using FastGaussQuadrature

# Independent scalar spectral functions from the Chebyshev expansion. Cut off
# before evaluating square roots; all scalar evaluation energies are in-band.
spectral_delta(n, x; delta=1e-5) = abs(x) < 1-delta ?
    2cos(n*acos(x))/(pi*sqrt(1-x^2)) : 0.0
spectral_R(n, x; delta=1e-5) = abs(x) < 1-delta ?
    -2im*exp(-im*n*acos(x))/sqrt(1-x^2) : 0.0im
spectral_A(n, x; delta=1e-5) = conj(spectral_R(n,x; delta=delta))

@testset "Green functions and differential optical response" begin
    for x in (-0.4, 0.0, 0.3), n in 0:3
        @test KPMru.gn_R(x,n) ≈ spectral_R(n,x)
        @test KPMru.gn_A(x,n) ≈ conj(KPMru.gn_R(x,n))
        @test KPMru.Δn(x,n) ≈ spectral_delta(n,x)
        @test KPMru.gn_R(x,n)-KPMru.gn_A(x,n) ≈ -2pi*im*KPMru.Δn(x,n)
    end
    NC = 4
    v = ComplexF64[1, 0.2im, 0.3, -0.1im]
    M = reshape(ComplexF64.(1:NC^2), NC,NC) ./ NC^2
    w = 0.12
    g = [(n==0 ? 1 : 2)*KPMru.JacksonKernel(n,NC) for n in 0:NC-1]
    for x in (-0.3, 0.1)
        ref1 = sum(g[n+1]*v[n+1]*spectral_delta(n,x) for n in 0:NC-1)
        ref2 = sum(g[n+1]*g[m+1]*M[n+1,m+1] *
                   (spectral_R(n,x+w)*spectral_delta(m,x) +
                    spectral_delta(n,x)*spectral_A(m,x-w))
                   for n in 0:NC-1, m in 0:NC-1)
        @test KPMru.d_optical_cond1(v,NC,x) ≈ ref1
        @test KPMru.d_optical_cond2(M,NC,w,x) ≈ ref2
    end
    for (grid, values, scalar) in (
        (KPMru.d_optical_cond1(v,NC; e_range=[-0.5,0.5])...,
         x -> KPMru.d_optical_cond1(v,NC,x)),
        (KPMru.d_optical_cond2(M,NC,w; e_range=[-0.5,0.5])...,
         x -> KPMru.d_optical_cond2(M,NC,w,x)))
        @test length(grid) == 2NC
        @test values ≈ scalar.(grid)
        @test all(isfinite,values)
    end

    # Integrated optical APIs must run too, not just the differential variants.
    nodes, weights = gausschebyshevt(8NC)
    occupation = [x < 0 ? 1.0 : 0.0 for x in nodes]
    ref1 = -im/w * sum(weights .* occupation .*
                      [KPMru.d_optical_cond1(v,NC,x) for x in nodes])
    ref2 = -im/w * sum(weights .* occupation .*
                      [KPMru.d_optical_cond2(M,NC,w,x) for x in nodes])
    @test KPMru.optical_cond1(v,NC,w) ≈ ref1
    @test KPMru.optical_cond2(M,NC,w) ≈ ref2
end

@testset "CPGE nonzero contraction and integration" begin
    NC = 4
    w, Omega, delta = 0.12, 0.006, 0.02
    w1, w2 = w, Omega-w
    # A single nonzero coefficient isolates a tensor component and prevents a
    # zero-valued output from satisfying all the regression checks accidentally.
    Gamma = zeros(ComplexF64,NC,NC,NC)
    Gamma[2,3,1] = 0.8+0.3im
    before = copy(Gamma)
    n,m,p = 1,2,0
    h(j) = j==0 ? 1 : 2
    kernel_weight = prod(h(j)*KPMru.JacksonKernel(j,NC) for j in (n,m,p))
    integrand(x) = spectral_R(n,x+w1+w2; delta=delta)*spectral_R(m,x+w2; delta=delta)*spectral_delta(p,x; delta=delta) +
                   spectral_R(n,x+w1; delta=delta)*spectral_delta(m,x; delta=delta)*spectral_A(p,x-w2; delta=delta) +
                   spectral_delta(n,x; delta=delta)*spectral_A(m,x-w1; delta=delta)*spectral_A(p,x-w1-w2; delta=delta)
    for x in (-0.3,0.1)
        reference = Gamma[2,3,1]*kernel_weight*integrand(x)
        @test KPMru.d_cpge(Gamma,NC,w1,w2,x; δ=delta) ≈ reference
    end
    nodes, weights = gausschebyshevt(8NC)
    integral = sum(weights[i]*integrand(x) for (i,x) in enumerate(nodes) if x < 0)
    reference = Gamma[2,3,1]*kernel_weight*integral*im*Omega/(w1*w2)
    result = KPMru.cpge(Gamma,NC,w; Ω=Omega, δ=delta)
    @test isfinite(result)
    @test abs(reference) > 1e-6
    @test result ≈ reference atol=1e-10
    @test KPMru.cpge(2Gamma,NC,w; Ω=Omega, δ=delta) ≈ 2result
    @test Gamma == before
    @test_throws AssertionError KPMru.cpge(Gamma,NC,w; Ω=w)
    grid, values = KPMru.d_cpge(Gamma,NC,w1,w2; e_range=[-0.4,0.4], δ=delta)
    @test length(grid) == 2NC
    @test values ≈ KPMru.d_cpge.([Gamma],NC,w1,w2,grid; δ=delta)
    # Exercise the adaptive integration API on a smooth, cutoff-free subinterval.
    quad(f) = KPMru.quadgk(f,-0.4,-0.1)
    @test KPMru.Λnmp((n,m,p),w1,w2; δ=delta,quad=quad) ≈ first(quad(integrand))
    @test KPMru.Λn((0,); quad=quad) ≈
          2/pi*(asin(-0.1)-asin(-0.4))
end

@testset "DC response compared with spectral matrix sums" begin
    (; H,Ja,psi,T,NH,NC) = small_system()
    a, energy = 1.5, 0.2
    mu = KPMru.kpm_2d(H,Ja,Ja,NC,NH,NH; psi_in=psi)
    for E in (0.0,energy)
        coeff = [(n==0 ? 1 : 2)*KPMru.JacksonKernel(n,NC)*cos(n*acos(E/a)) for n in 0:NC-1]
        F = sum(coeff[n]*T[n] for n in 1:NC)
        reference = tr(F*Ja*F*Ja)/(NH*a)
        @test KPMru.dc_cond_single(mu,a,E) ≈ reference atol=1e-12
        if E==0
            @test KPMru.dc_cond0(mu,a) ≈ reference atol=1e-12
        end
    end
    orders = [4,8,6]
    for Ef in (0.0,energy)
        reference = map(orders) do order
            F = sum((n==0 ? 1 : 2)*KPMru.JacksonKernel(n,order)*
                    cos(n*acos(Ef/a))*T[n+1] for n in 0:order-1)
            tr(F*Ja*F*Ja)/(NH*a)
        end
        @test KPMru.dc_long(H,Ja,a,orders,NH,NH; psi_in=psi,Ef=Ef,debug_mode=false) ≈ reference
    end
    all_samples = KPMru.dc_long(H,Ja,a,orders,NH,NH; psi_in=psi,avg_NR=false,debug_mode=false)
    average = KPMru.dc_long(H,Ja,a,orders,NH,NH; psi_in=psi,debug_mode=false)
    @test vec(sum(all_samples; dims=2)) ./ NH ≈ average
    h = 1e-5
    derivative = (KPMru.d_dc_cond(mu,a,energy+h)-KPMru.d_dc_cond(mu,a,energy-h))/(2h)
    @test KPMru.d_dc_cond(mu,a,energy; dE_order=1) ≈ derivative rtol=1e-7
    @test KPMru.d_dc_cond(mu,a,[-2.0,2.0]) == zeros(2)
end
