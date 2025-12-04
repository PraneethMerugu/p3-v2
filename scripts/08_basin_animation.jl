using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using DynamicalSystems
using GLMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using Logging

# 1. Define the Stroboscopic Map Integration
# This function advances the state by exactly one period T
function duffing_step!(u_new, u, p, t)
    # p = [μ, ε, A, ω]
    ω = p[4]
    T = 2π / ω

    # Convert u to SVector
    u_svec = SVector{2}(u)

    # Callback to terminate if trajectory diverges
    condition(u, t, integrator) = norm(u) - 100.0
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(duffing_rule, u_svec, (0.0, T), p, callback=cb)

    # Suppress warnings and use relaxed tolerances for speed
    sol = with_logger(NullLogger()) do
        solve(prob, Rodas5P(), save_everystep=false, maxiters=1000,
            abstol=1e-4, reltol=1e-4, verbose=false)
    end

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Terminated
        u_new .= sol[end]
    else
        # If failed, assume divergent
        u_new .= [1e3, 1e3]
    end
    return nothing
end

# 2. Visualization & Computation Setup
function animate_fractal_basins(output_path)
    println("Initializing Basin Animation...")

    # Parameters
    # We will scan Amplitude A from 0.1 (Period 1) to 0.35 (Chaos/Multi-period)
    A_range = range(0.3, 0.5, length=20)

    # Fixed parameters [μ, ε, A_placeholder, ω]
    p_init = [-0.225, 0.25, A_range[1], 1.0]

    # Define the discrete system (Stroboscopic Map)
    # u0 is a dummy initial state
    ds_map = DeterministicIteratedMap(duffing_step!, [0.1, 0.1], p_init)

    # Define the Grid of Initial Conditions
    # Resolution: 100x100 is a good balance for animation speed vs detail
    xg = range(-2.5, 2.5, length=100)
    yg = range(-2.0, 2.0, length=100)
    grid = (xg, yg)

    # Initialize Mapper
    # 'AttractorsViaRecurrences' finds attractors automatically.
    # sparse=false returns a full matrix of IDs for the grid.
    mapper = AttractorsViaRecurrences(ds_map, grid; sparse=false)

    # 3. Setup Makie Plot
    fig = Figure(size=(800, 700))
    ax = Axis(fig[1, 1], title="Fractal Basins of Attraction", xlabel="x", ylabel="y")

    # Initial Calculation
    println("Computing initial frame...")
    basins, attractors = basins_of_attraction(mapper, grid)

    # Create an Observable for the heatmap data
    # basins is a Matrix{Int} where each Int is an attractor ID (1, 2, -1 for divergent, etc.)
    basins_obs = Observable(basins)

    # Custom colormap to handle discrete attractor IDs
    # -1 (divergent) is usually colored black or white
    heatmap!(ax, xg, yg, basins_obs, colormap=:glasbey_hv_n256)

    # 4. Animation Loop
    println("Starting render loop...")
    record(fig, output_path, A_range; framerate=5) do A_val
        # Update Parameter in the System
        # Index 3 is Amplitude A
        ds_map.p[3] = A_val

        # Update Title
        ax.title = "Forcing A = $(round(A_val, digits=3))"

        # IMPORTANT: Reset the mapper for the new system state!
        # This clears the cache of found attractors so we find new ones for the new parameter.
        # We assume attractors change, so we don't map previous attractors to new ones here (simplest approach).
        mapper = AttractorsViaRecurrences(ds_map, grid; sparse=false)

        # Compute Basins
        new_basins, new_attractors = basins_of_attraction(mapper, grid)

        # Update Plot
        basins_obs[] = new_basins
    end

    println("Animation saved to $output_path")
end

# Run
output_file = joinpath(@__DIR__, "../output/fractal_basins.mp4")
animate_fractal_basins(output_file)
