using Test

# The README's Julia fences form one CPU session. Installation and optional GPU
# commands are explicitly presented as text, outside this executable session.
readme = read(joinpath(@__DIR__, "..", "README.md"), String)
blocks = collect(eachmatch(r"(?mi)^ {0,3}```julia\h*\r?\n(.*?)^ {0,3}```\h*\r?$"s, readme))
@test length(blocks) >= 5
@test length(blocks) == count(r"(?mi)^ {0,3}```julia\h*\r?$", readme)

mktempdir() do sandbox_dir
    sandbox = Module(gensym(:ReadmeExamples))
    cd(sandbox_dir) do
        for (i, block) in enumerate(blocks)
            @testset "README Julia block $i" begin
                snippet = joinpath(sandbox_dir, "readme_$(i).jl")
                write(snippet, block.captures[1])
                @test begin
                    Base.include(sandbox, snippet)
                    true
                end
            end
        end
    end

    # Newly included module bindings participate in world age on newer Julia.
    # Read them in the latest world when checking the examples' useful outputs.
    Base.invokelatest() do
        @test length(sandbox.E) == length(sandbox.rho) == 51
        @test all(isfinite, sandbox.rho)
        @test sandbox.H_norm ≈ sandbox.H_orig / sandbox.a
        @test length(sandbox.mu) == sandbox.NC
        @test sandbox.mu[1] ≈ 1
        @test sandbox.E_grid ≈ collect(range(-1.5, 1.5; length=51))
        @test length(sandbox.rho_E) == 51
        @test all(isfinite, sandbox.rho_E)
        @test minimum(sandbox.rho_E) >= -1e-10
        @test isfinite(sandbox.rho_0)
        @test isfinite(sandbox.d2rho_0)
        @test size(sandbox.mu_2d_xx) == (sandbox.NC, sandbox.NC)
        @test length(sandbox.dsigma_xx) == length(sandbox.Evals)
        @test all(isfinite, sandbox.dsigma_xx)
        @test size(sandbox.mu_3d_xxx) == (sandbox.NC3, sandbox.NC3, sandbox.NC3)
        @test isfinite(sandbox.dchi_xxx)
    end
end
