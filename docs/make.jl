# Docs build script for Documenter.jl
# Render plots without opening GUI windows on CI or local docs builds.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using KPMru
using Documenter

makedocs(
    sitename = "KPMru",
    modules  = [KPMru],
    pages    = ["Home" => "index.md"],
    build    = joinpath("build", "dev"),
)
