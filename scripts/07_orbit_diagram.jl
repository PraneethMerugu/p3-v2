using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using DynamicalSystems
using GLMakie
using JLD2

# Explicitly include visualizations if needed
include(joinpath(@__DIR__, "../src/visualizations.jl"))
using .Visualizations

println("--- Orbit Diagram Analysis (Brute Force) ---")

# 1. Configuration
# Scan Amplitude A to see the route to chaos
p_scan = range(0.01, 0.6, length=400)
fixed_p = [-0.225, 0.25, 0.0, 1.0] # A is placeholder at index 3

# 2. Compute Loop
println("Computing Orbit Diagram over $(length(p_scan)) steps...")

output_data = Vector{Point2f}()
sizehint!(output_data, length(p_scan) * 100)

# Initialize System
ds = get_ds([0.1, 0.1], fixed_p)

# Identifying Indices
idx_A = 3
idx_ω = 4

println("Type of ds: ", typeof(ds))
println("Parameters p: ", ds.integ.p)
println("Index A: ", idx_A)
println("Index ω: ", idx_ω)


for A_val in p_scan
    # Update Parameter A
    set_parameter!(ds, idx_A, A_val)

    # Calculate Period T
    ω = ds.integ.p[idx_ω]
    T = 2π / ω

    # Burn-in (Transient removal)
    step!(ds.integ, 500 * T, true)

    # Stroboscopic Sampling
    # We collect 50 points per parameter value
    for _ in 1:50
        step!(ds.integ, T, true)
        # Store (A, x)
        push!(output_data, Point2f(A_val, ds.integ.u[1]))
    end
end

# 3. Plot
println("Rendering...")
fig = Figure(size=(1000, 600))
ax = Axis(fig[1, 1],
    xlabel="Forcing Amplitude (A)",
    ylabel="x (Stroboscopic)",
    title="Bifurcation Diagram (Brute Force)",
    limits=(nothing, nothing, -3, 3)
)

# Use high transparency to reveal density
scatter!(ax, output_data, markersize=1.0, color=(:black, 0.2))

save(joinpath(@__DIR__, "../output/orbit_diagram.png"), fig)
println("Saved orbit_diagram.png")
