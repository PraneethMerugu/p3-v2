using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using P3_V2
using DynamicalSystems
using Statistics

println("--- Debugging Data Ranges ---")

# 1. Check Stroboscopic Map (06)
println("\n1. Checking 06_stroboscopic_map.jl parameters...")
p_chaos = [-0.225, 0.25, 0.32, 1.0]
ds = get_ds([0.1, 0.1], p_chaos)
ω = p_chaos[4]
T = 2π / ω
println("Running trajectory for 06...")
try
    tr, t = trajectory(ds, 5000 * T; Δt=T, Ttr=1000 * T)
    x_min, x_max = extrema(tr[:, 1])
    y_min, y_max = extrema(tr[:, 2])
    println("06 Data Range:")
    println("  x: [$x_min, $x_max]")
    println("  y: [$y_min, $y_max]")
    if any(isnan, tr) || any(isinf, tr)
        println("  WARNING: NaNs or Infs detected in 06 data!")
    end
catch e
    println("  Error running 06 trajectory: ", e)
end

# 2. Check Orbit Diagram (07)
println("\n2. Checking 07_orbit_diagram.jl logic...")
# We'll simulate a smaller scan to save time, but cover the range
p_scan = range(0.2, 0.6, length=50)
fixed_p = [-0.225, 0.25, 0.0, 1.0]
ds = get_ds([0.1, 0.1], fixed_p)
idx_A = 3
idx_ω = 4
output_data_x = Float64[]
output_data_A = Float64[]

println("Running scan for 07...")
for A_val in p_scan
    set_parameter!(ds, idx_A, A_val)
    ω = ds.integ.p[idx_ω]
    T = 2π / ω

    try
        step!(ds.integ, 100 * T, true) # Shorter burn-in
        for _ in 1:10
            step!(ds.integ, T, true)
            push!(output_data_x, ds.integ.u[1])
            push!(output_data_A, A_val)
        end
    catch e
        println("  Error at A=$A_val: ", e)
    end
end

if !isempty(output_data_x)
    x_min, x_max = extrema(output_data_x)
    println("07 Data Range (x):")
    println("  x: [$x_min, $x_max]")
else
    println("07 produced no data.")
end
