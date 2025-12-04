using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using GLMakie, ModelingToolkit, DynamicalSystems

include(joinpath(@__DIR__, "../src/visualizations.jl"))
using .Visualizations

Visualizations.set_project_theme!()

println("Starting Interactive App...")

# 1. Setup System
p_syms = parameters(duffing_sys)
val_map = Dict(μ => -0.225, ε => 0.25, A => 0.3, ω => 1.0)
p_vals = [val_map[p] for p in p_syms]
u0 = [0.1, 0.1]

# Create CoupledODEs
ds = get_coupled_odes(duffing_sys, u0, p_vals)

# Find indices
idx_ε = ModelingToolkit.parameter_index(duffing_sys, ε)
idx_A = ModelingToolkit.parameter_index(duffing_sys, A)
idx_ω = ModelingToolkit.parameter_index(duffing_sys, ω)

# 2. Visualization Setup
fig = Figure(size=(1000, 800))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="Duffing Oscillator Explorer")
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
    # Update parameters using ParameterIndex
    ds.integ.p[idx_ε] = e
    ds.integ.p[idx_A] = a
    ds.integ.p[idx_ω] = w

    # Re-initialize and run
    # We use reinit! to reset state, or step! to continue?
    # For exploring attractors, reinit! is safer to see the attractor for *this* parameter set.
    # But for "tracking", step! is better.
    # Let's use reinit! with a fixed random start or keep previous?
    # Let's keep previous state to show hysteresis!
    # But if it diverges, we might need reset.

    # For this visualizer, let's reinit! to u0 to be deterministic.
    reinit!(ds.integ, u0)

    # Run trajectory
    tr, _ = trajectory(ds, 200.0; Ttr=500.0) # Transient is important
    return Point2f.(tr[:, 1], tr[:, 2])
end

lines!(ax, points, color=:cyan, linewidth=1.5)

display(fig)
println("App Running.")
