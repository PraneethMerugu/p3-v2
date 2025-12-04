using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using GLMakie, JLD2, P3_V2

# Visualizations are now exported by P3_V2

println("--- Chaos & Complexity Analysis (Plot) ---")

input_file = joinpath(@__DIR__, "../output/chaos_data.jld2")
if !isfile(input_file)
    error("Data file not found: $input_file. Run 03_chaos_compute.jl first.")
end

println("Loading data from $input_file...")
data = load(input_file)
lz_matrix = data["lz_matrix"]
A_range = data["A_range"]
ω_range = data["ω_range"]

# Plot
println("Plotting Heatmap...")
plot_heatmap(A_range, ω_range, lz_matrix, "Lempel-Ziv Complexity", joinpath(@__DIR__, "../output/chaos_heatmap.png"))
