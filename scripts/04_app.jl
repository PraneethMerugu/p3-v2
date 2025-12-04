using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using GLMakie, DynamicalSystems

# Visualizations are now exported by P3_V2

set_project_theme!()

println("Starting Interactive App...")

# 1. Setup System
# p = [μ, ε, A, ω]
# μ = -0.225, ε = 0.25, A = 0.3, ω = 1.0
p_vals = [-0.225, 0.25, 0.3, 1.0]
u0 = [0.1, 0.1]

# Create CoupledODEs
ds = get_ds(u0, p_vals)

# Parameter indices
# p = [μ, ε, A, ω]
idx_ε = 2
idx_A = 3
idx_ω = 4

# 2. Visualization Setup
fig = Figure(size=(1000, 800))
ax = Axis(fig[1, 1], title="Duffing Oscillator Phase Space", xlabel="x", ylabel="y")
xlims!(ax, -2.5, 2.5)
ylims!(ax, -2.0, 2.0)

# Sliders
sg = SliderGrid(fig[2, 1],
    (label="Damping (ε)", range=0:0.001:0.5, format="{:.3f}", startvalue=0.25),
    (label="Forcing Amp (A)", range=0:0.01:1.5, format="{:.2f}", startvalue=0.3),
    (label="Frequency (ω)", range=0.5:0.01:2.0, format="{:.2f}", startvalue=1.0)
)

ε_obs = sg.sliders[1].value
A_obs = sg.sliders[2].value
ω_obs = sg.sliders[3].value

# 3. Reactive Trajectory
points = lift(ε_obs, A_obs, ω_obs) do e, a, w
    # Update parameters using integer indices
    set_parameter!(ds, idx_ε, e)
    set_parameter!(ds, idx_A, a)
    set_parameter!(ds, idx_ω, w)

    # Re-initialize and run
    reinit!(ds, u0)

    # Run trajectory
    tr, _ = trajectory(ds, 200.0; Ttr=500.0) # Transient is important
    return Point2f.(tr[:, 1], tr[:, 2])
end

# Add dynamic vector field
# We pass the slider observables directly
plot_vector_field!(ax, -0.225, ε_obs, A_obs, ω_obs)

lines!(ax, points, color=:cyan, linewidth=1.5)

display(fig)
println("App Running.")
