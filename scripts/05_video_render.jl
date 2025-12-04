using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using GLMakie, ModelingToolkit, Attractors

include(joinpath(@__DIR__, "../src/visualizations.jl"))
using .Visualizations

set_project_theme!()

println("--- Video Rendering Studio ---")

# 1. Animation B: Limit Cycle Collision
render_collision() = animate_collision(duffing_sys, range(0.23, 0.20, length=100), joinpath(@__DIR__, "../output/homoclinic_collision.mp4"))

# 2. Animation D: Sensitivity (Butterfly Effect)
render_sensitivity() = animate_sensitivity(duffing_sys, joinpath(@__DIR__, "../output/sensitivity.mp4"))

# Run
render_collision()
render_sensitivity()
