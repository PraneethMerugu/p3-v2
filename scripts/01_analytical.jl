using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2

# Define parameters
# μ = -0.225, ε = 0.225 (Hopf), A = 0, ω = 1.0
# p = [μ, ε, A, ω]
p_val = [-0.225, 0.225, 0.0, 1.0]

# Equilibrium points
# u = [x, y]
P0 = [0.0, 0.0]
P1 = [1.0, 0.0]
P2 = [-1.0, 0.0]

println("--- Stability Analysis (ε = 0.225) ---")

# Analyze P0
res0 = analyze_equilibrium(P0, p_val)
println("\nP0 (0,0): ", res0.classification)
println("Eigenvalues: ", res0.eigenvalues)
println("Jacobian:\n", res0.jacobian)

# Analyze P1
res1 = analyze_equilibrium(P1, p_val)
println("\nP1 (1,0): ", res1.classification)
println("Eigenvalues: ", res1.eigenvalues)

# Analyze P2
res2 = analyze_equilibrium(P2, p_val)
println("\nP2 (-1,0): ", res2.classification)
println("Eigenvalues: ", res2.eigenvalues)

println("\n--- Analysis Complete ---")
