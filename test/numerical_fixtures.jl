using LinearAlgebra, Random, SparseArrays

# A full basis makes the trace exact, removing stochastic tolerance/flakiness.
# Keep its columns uncentered and unit-normalized (do not use force_norm).
function small_system(; NC=8)
    rng = MersenneTwister(405)
    hermitian() = begin
        A = randn(rng, ComplexF64, 6, 6)
        (A + A') / 2
    end
    H = hermitian()
    H .*= 0.7 / opnorm(H)
    Ja, Jb, Jc = hermitian(), hermitian(), hermitian()
    psi = Matrix{ComplexF64}(I, 6, 6)
    eig = eigen(Hermitian(H))
    T = [eig.vectors * Diagonal(cos.(n .* acos.(eig.values))) * eig.vectors'
         for n in 0:NC-1]
    return (; H, Ja, Jb, Jc, psi, T, NH=6, NC)
end
