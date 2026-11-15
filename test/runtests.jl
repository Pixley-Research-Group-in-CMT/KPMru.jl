using KPMru
using Test

include("numerical_fixtures.jl")

@testset "Numerical APIs" begin
    include("moments_test.jl")
    include("dos_test.jl")
    include("utilities_test.jl")
    include("response_test.jl")
end

@testset "README examples" begin
    include("readme_test.jl")
end

@testset "kernels/kernels.jl" begin
    include("test_kernel.jl")
end

@testset "util/chebyshev_iteration.jl" begin
    include("chebyshev_iteration_test.jl")
end

@testset "KPMru.jl" begin
    include("integration_test.jl")
end

@testset "optional CUDA extension migration" begin
    include("optional_cuda_test.jl")
end
