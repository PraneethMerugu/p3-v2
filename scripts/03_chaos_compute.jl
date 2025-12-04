using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using JLD2
using DynamicalSystems

println("--- Chaos & Complexity Analysis (Compute) ---")

# 1. Setup System
# Chaos parameters: μ = -0.225, ε = 0.25, A = 0.3, ω = 1.0
# p = [μ, ε, A, ω]
p_chaos = [-0.225, 0.25, 0.3, 1.0]
u0 = [0.1, 0.1]

ds = get_ds(u0, p_chaos)

# 2. Compute Metrics
println("Computing Lempel-Ziv Complexity...")
lz = compute_lz_complexity(ds, 5000.0; Ttr=1000.0)
println("LZ Complexity (A=0.3): ", lz)

println("Computing Markov Matrix...")
M = compute_markov_matrix(ds, 5000.0; Ttr=1000.0)
println("Markov Matrix (First 5 rows):\n", M[1:min(5, size(M, 1)), :])

# 3. Heatmap Scan
println("Computing Chaos Heatmap (Grid Search)...")

# Define grid
A_range = range(0.01, 1.51, length=200)
ω_range = range(0.01, 2.01, length=200)

# Fixed parameters
fixed_p = Dict(:μ => -0.225, :ε => 0.25)

# Compute
# This uses the new parallelized function in src/chaos.jl
# We pass nothing for duffing_sys as it's no longer used
lz_matrix = compute_chaos_grid(nothing, A_range, ω_range, fixed_p)

# Save results
output_file = joinpath(@__DIR__, "../output/chaos_data.jld2")
save(output_file, "lz_matrix", lz_matrix, "A_range", A_range, "ω_range", ω_range)
println("Saved computation results to $output_file")
