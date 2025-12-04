# Refactoring Plan: Duffing Oscillator Analysis via DynamicalSystems.jl

**Objective:**
Refactor the codebase to exclusively use `DynamicalSystems.jl` for simulation and analysis, removing all dependencies on `BifurcationKit.jl`. This aligns with the project goal of analyzing chaos (stroboscopic maps, complexity) rather than tracking unstable equilibria via continuation.

**Difficulty:** Low to Medium.
**Key Change:** Replace symbolic `ModelingToolkit` definitions with standard Julia functions using `StaticArrays` for zero-allocation performance.

---

## 1. Refactor `src/model.jl`
**Action:** Replace the symbolic variables with a high-performance function definition.

```julia
module Model
    using StaticArrays

    # Export parameters for consistency across scripts
    export duffing_rule, p_default

    # Default parameter set
    # p = [μ, ε, A, ω]
    const p_default = [
        -0.225, # μ (Nonlinear Damping)
        0.25,   # ε (Linear Damping)
        0.3,    # A (Forcing Amplitude)
        1.0     # ω (Forcing Frequency)
    ]

    """
    duffing_rule(u, p, t)

    The equations of motion:
    x' = y
    y' = x - x^3 - ε*y - μ*x^2*y + A*cos(ω*t)
    """
    function duffing_rule(u, p, t)
        x, y = u
        μ, ε, A, ω = p
        
        dx = y
        dy = x - x^3 - ε*y - μ*x^2*y + A*cos(ω*t)
        
        return SVector(dx, dy)
    end
end
```

---

## 2. Refactor `src/chaos.jl`
**Action:** Simplify the system initialization to use `CoupledODEs` directly.

```julia
module Chaos
    using DynamicalSystems
    using ..Model # Import the rule from above

    export get_ds

    """
    get_ds(u0, p)

    Returns a CoupledODEs (ContinuousDynamicalSystem) initialized with the Duffing rule.
    """
    function get_ds(u0=[0.1, 0.1], p=p_default)
        # CoupledODEs is the standard struct for continuous systems in DynamicalSystems.jl v3+
        return CoupledODEs(duffing_rule, u0, p)
    end
end
```

---

## 3. Create `scripts/07_orbit_diagram.jl`
**Action:** Implement the "Brute Force" bifurcation diagram. This replaces `BifurcationKit`'s continuation method.

```julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
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
p_scan = range(0.2, 0.6, length=400) 
fixed_p = [-0.225, 0.25, 0.0, 1.0] # A is placeholder at index 3

# 2. Compute Loop
println("Computing Orbit Diagram over $(length(p_scan)) steps...")

output_data = Vector{Point2f}()
sizehint!(output_data, length(p_scan) * 100)

# Initialize System
ds = Chaos.get_ds([0.1, 0.1], fixed_p)

# Identifying Indices
idx_A = 3
idx_ω = 4

for A_val in p_scan
    # Update Parameter A
    set_parameter!(ds, idx_A, A_val)
    
    # Calculate Period T
    ω = current_parameter(ds)[idx_ω]
    T = 2π / ω
    
    # Burn-in (Transient removal)
    step!(ds, 500 * T)
    
    # Stroboscopic Sampling
    # We collect 50 points per parameter value
    for _ in 1:50
        step!(ds, T)
        # Store (A, x)
        push!(output_data, Point2f(A_val, current_state(ds)[1]))
    end
end

# 3. Plot
println("Rendering...")
fig = Figure(size=(1000, 600))
ax = Axis(fig[1,1], 
    xlabel = "Forcing Amplitude (A)", 
    ylabel = "x (Stroboscopic)", 
    title = "Bifurcation Diagram (Brute Force)"
)

# Use high transparency to reveal density
scatter!(ax, output_data, markersize = 1.0, color = (:black, 0.2))

save(joinpath(@__DIR__, "../output/orbit_diagram.png"), fig)
println("Saved orbit_diagram.png")
```

---

## 4. Refactor `scripts/06_stroboscopic_map.jl`
**Action:** Fix the `Delta_t` keyword issue by using the correct `Δt`.

```julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using P3_V2
using DynamicalSystems
using GLMakie

println("--- Generating Stroboscopic Map ---")

# 1. Setup (Chaotic Regime A=0.32)
p_chaos = [-0.225, 0.25, 0.32, 1.0] 
ds = Chaos.get_ds([0.1, 0.1], p_chaos)

# 2. Define Sampling Period
ω = p_chaos[4]
T = 2π / ω
println("Sampling period T = $T")

# 3. Generate Orbit
# Δt = T ensures we sample exactly once per period (Stroboscopic)
tr, t = trajectory(ds, 5000 * T; Δt = T, Ttr = 1000 * T)

# 4. Visualization
fig = Figure(size=(800,600))
ax = Axis(fig[1,1], 
    title="Stroboscopic Map (A=0.32)", 
    xlabel="x", ylabel="y"
)

# Plot the points
scatter!(ax, tr[:, 1], tr[:, 2], markersize=2, color=:red, label="Stroboscopic Points")

# Add background flow (faint)
tr_flow, _ = trajectory(ds, 20 * T; Δt = 0.01, Ttr = 0)
lines!(ax, tr_flow[:, 1], tr_flow[:, 2], color=(:grey, 0.3), linewidth=0.5)

axislegend(ax)
save(joinpath(@__DIR__, "../output/stroboscopic_map.png"), fig)
println("Map saved.")
```

---

## 5. Cleanup
**Action:** Delete `scripts/02_bifurcation.jl` and any other scripts referencing `BifurcationKit` to prevent confusion.