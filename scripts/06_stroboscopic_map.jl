using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using DynamicalSystems
using GLMakie

println("--- Generating Stroboscopic Map ---")

# 1. Setup (Chaotic Regime A=0.32)
p_chaos = [-0.225, 0.25, 0.02, 1.0]
ds = get_ds([0.1, 0.1], p_chaos)

# 2. Define Sampling Period
ω = p_chaos[4]
T = 2π / ω
println("Sampling period T = $T")

# 3. Generate Orbit
# Δt = T ensures we sample exactly once per period (Stroboscopic)
tr, t = trajectory(ds, 5000 * T; Δt=T, Ttr=1000 * T)

# 4. Visualization
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1],
    title="Stroboscopic Map (A=0.32)",
    xlabel="x", ylabel="y",
    limits=(-3, 3, -3, 3)
)

# Plot the points
scatter!(ax, tr[:, 1], tr[:, 2], markersize=2, color=:red, label="Stroboscopic Points")

# Add background flow (faint)
tr_flow, _ = trajectory(ds, 20 * T; Δt=0.01, Ttr=0)
lines!(ax, tr_flow[:, 1], tr_flow[:, 2], color=(:grey, 0.3), linewidth=0.5)

axislegend(ax)
save(joinpath(@__DIR__, "../output/stroboscopic_map.png"), fig)
println("Map saved.")
