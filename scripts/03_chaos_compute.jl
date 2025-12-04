using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using ModelingToolkit, JLD2

println("--- Chaos & Complexity Analysis (Compute) ---")

# 1. Setup System
p_syms = parameters(duffing_sys)
# Chaos parameters: μ = -0.225, ε = 0.25, A = 0.3, ω = 1.0
val_map = Dict(μ => -0.225, ε => 0.25, A => 0.3, ω => 1.0)
p_vals = [val_map[p] for p in p_syms]

u0 = [0.1, 0.1]
ds = get_coupled_odes(duffing_sys, u0, p_vals)

# 2. Compute Metrics
println("Computing Lempel-Ziv Complexity...")
lz = compute_lz_complexity(ds, 5000.0; Ttr=1000.0)
println("LZ Complexity (A=0.3): ", lz)

println("Computing Markov Matrix...")
M = compute_markov_matrix(ds, 5000.0; Ttr=1000.0)
println("Markov Matrix:\n", M)

# 3. Heatmap Scan
println("Computing Chaos Heatmap (Grid Search)...")

# Define grid
A_range = range(0.0, 1.5, length=50)
ω_range = range(0.5, 2.0, length=50)

# Fixed parameters
fixed_p = Dict(μ => -0.225, ε => 0.25)

# Compute
# This uses the new parallelized function in src/chaos.jl
lz_matrix = compute_chaos_grid(duffing_sys, A_range, ω_range, fixed_p)

# Save results
output_file = joinpath(@__DIR__, "../output/chaos_data.jld2")
save(output_file, "lz_matrix", lz_matrix, "A_range", A_range, "ω_range", ω_range)
println("Saved computation results to $output_file")
