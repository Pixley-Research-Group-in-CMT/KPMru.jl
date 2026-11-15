using Statistics

@testset "Scaling, vector operations and scalar utilities" begin
    s = small_system()
    original = sparse(s.H)
    a, Hn = KPMru.normalizeH(original; fixed_a=2.0)
    @test a == 2.0
    @test Hn ≈ original/2
    @test original == sparse(s.H)
    @test KPMru.wrapAdd(1, -1, 5) == 5
    @test KPMru.wrapAdd(5, 1, 5) == 1
    @test KPMru.isNotBoundary([1,2], [4,2], [4,4], [1,0]) == 0
    @test KPMru.isNotBoundary([1,2], [4,2], [4,4], [0,0]) == 1
    @test KPMru.fermiFunction(0.2, 0.2, 10.0) == 0.5
    @test KPMru.fermiFunction(; E=0, E_f=0, beta=3) == 0.5
    for beta in (2.0, Inf)
        f = KPMru.fermiFunctions(0.0, beta)
        @test f(0.0) == 0.5
        @test f(-0.3) ≈ 1-f(0.3)
        @test 0 <= f(0.5) <= f(-0.5) <= 1
    end
    @test KPMru.fermiFunctions(; E_f=0, beta=Inf).([-1.0,0.0,1.0]) == [1.0,0.5,0.0]
    V = ComplexF64[1 2; 2 -1; 4 3]
    centered = copy(V)
    KPMru.normalize_by_col(centered, 2)
    @test vec(sum(centered; dims=1)) ≈ zeros(2) atol=1e-15
    @test [norm(centered[:,i]) for i in 1:2] ≈ ones(2)
    uncentered = copy(V)
    KPMru.normalize_by_col(uncentered, 2; centering=false)
    @test uncentered ≈ hcat([V[:,i]/norm(V[:,i]) for i in 1:2]...)
    Q = KPMru.gram_schmidt(V)
    @test Q'Q ≈ Matrix{ComplexF64}(I,2,2) atol=1e-14
    @test V == ComplexF64[1 2; 2 -1; 4 3]
    @test Q*Q'*V ≈ V
    KPMru.gram_schmidt!(V)
    @test V ≈ Q
end

@testset "Chebyshev transforms and kernel weights" begin
    n = collect(0:7)
    x = [-0.8, 0.0, 0.3]
    coeff = ComplexF64[1, 2im, 0.3, -0.2im, 0.4, 0.1, 0.0, -0.3]
    direct = [sum(coeff[j+1]*cos(j*acos(t)) for j in n) for t in x]
    @test KPMru.chebyshev_lin_trans(x,n,coeff) ≈ direct
    @test KPMru.chebyshev_lin_trans(x[1],n,coeff) ≈ direct[1]
    for j in n, t in x
        @test KPMru.chebyshevT_accurate(j,t) ≈ cos(j*acos(t)) atol=1e-14
    end
    NC = length(n)
    mu = reshape(ComplexF64.(1:NC^2),NC,NC)
    before = copy(mu)
    kernel = KPMru.JacksonKernel
    weights = [(j==0 ? 1 : 2)*kernel(j,NC) for j in n]
    for dims in (Int[], [1], [2], [1,2])
        expected = [mu[i,j] * (1 in dims ? weights[i] : 1) *
                    (2 in dims ? weights[j] : 1) for i in 1:NC, j in 1:NC]
        @test KPMru.mu2D_apply_kernel_and_h(mu,NC,kernel; dims=dims) ≈ expected
        @test KPMru.muND_apply_kernel_and_h_no_mutate(mu,NC,kernel; dims=dims) ≈ expected
    end
    @test mu == before
    @test KPMru.mu2D_apply_kernel_and_h_dims1(mu,NC,weights) ≈ mu .* weights
    @test KPMru.mu2D_apply_kernel_and_h_dims2(mu,NC,weights) ≈ mu .* weights'
end
