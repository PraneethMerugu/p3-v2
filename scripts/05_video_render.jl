using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using GLMakie

# Visualizations are now exported by P3_V2

set_project_theme!()

println("--- Video Rendering Studio ---")

# 1. Animation B: Limit Cycle Collision
# animate_collision(ε_range, output_path)
render_collision() = animate_collision(range(0.23, 0.20, length=100), joinpath(@__DIR__, "../output/homoclinic_collision.mp4"))

# 2. Animation D: Sensitivity (Butterfly Effect)
# animate_sensitivity(output_path)
render_sensitivity() = animate_sensitivity(joinpath(@__DIR__, "../output/sensitivity.mp4"))

# Run
render_collision()
render_sensitivity()
