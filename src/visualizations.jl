module Visualizations

using GLMakie, DynamicalSystems, DifferentialEquations, LinearAlgebra
using P3_V2
using P3_V2.Model
using P3_V2.Chaos

export animate_collision, animate_sensitivity, plot_heatmap, set_project_theme!

function set_project_theme!()
    set_theme!(theme_black())
    update_theme!(
        fontsize=20,
        Axis=(
            xgridvisible=false,
            ygridvisible=false,
        )
    )
end

"""
    plot_nullclines!(ax, x_range, μ, ε, A, ω, t; linewidth=2.0)

Plots the x and y nullclines for the Modified Duffing Oscillator.
Handles singularities in the y-nullcline where the denominator (ε + μ*x^2) is near zero.
"""
function plot_nullclines!(ax, x_range, μ, ε, A, ω, t; linewidth=2.0)
    # 1. x-nullcline: y = 0
    hlines!(ax, [0.0], color=:orange, linestyle=:dash, linewidth=linewidth, label="x-nullcline")

    # 2. y-nullcline: y = (x - x^3 + A*cos(ω*t)) / (ε + μ*x^2)
    y_vals = Float64[]
    valid_x = Float64[]

    for x_val in x_range
        denom = ε + μ * x_val^2

        # Avoid division by zero or extreme values near asymptotes
        if abs(denom) > 1e-4
            num = x_val - x_val^3 + A * cos(ω * t)
            y_val = num / denom

            # Clamp for visualization purposes if needed, or just let it go off-screen
            # Here we just check if it's reasonable to plot to avoid artifacts
            if abs(y_val) < 100.0
                push!(y_vals, y_val)
                push!(valid_x, x_val)
            else
                push!(y_vals, NaN)
                push!(valid_x, x_val)
            end
        else
            push!(y_vals, NaN)
            push!(valid_x, x_val)
        end
    end

    lines!(ax, valid_x, y_vals, color=:cyan, linestyle=:dash, linewidth=linewidth, label="y-nullcline")
end

"""
    animate_collision(sys, ε_range, output_path)

Generates an animation showing the limit cycle collision as ε varies.
"""
function animate_collision(ε_range, output_path)
    println("Rendering Collision Animation to $output_path...")

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], title="Homoclinic Collision", xlabel="x", ylabel="y")
    xlims!(ax, -1.5, 1.5)
    ylims!(ax, -1.0, 1.0)

    # Fixed parameters
    μ_val = -0.225
    A_val = 0.0
    ω_val = 1.0

    record(fig, output_path, ε_range; framerate=10) do e_val
        ax.title = "ε = $(round(e_val, digits=4))"
        empty!(ax)

        # Plot Nullclines
        plot_nullclines!(ax, range(-1.5, 1.5, length=200), μ_val, e_val, A_val, ω_val, 0.0)

        # Define parameters for this frame
        # p = [μ, ε, A, ω]
        p_vals = [μ_val, e_val, A_val, ω_val]

        # Simulate from P1 (1.0, 0.0) to find the attractor
        u0 = [1.0, 0.0]
        ds = get_ds(u0, p_vals)

        # Trajectory
        tr, _ = trajectory(ds, 200.0; Ttr=500.0)
        lines!(ax, tr[:, 1], tr[:, 2], color=:green, linewidth=2, label="Trajectory")

        # Plot Saddle P0 (0,0)
        scatter!(ax, [0], [0], color=:red, markersize=15, label="Saddle (0,0)")

        # Add Legend
        axislegend(ax, position=:rt)
    end
    println("Saved $output_path")
end

"""
    animate_sensitivity(output_path)

Generates an animation demonstrating the butterfly effect (sensitivity to initial conditions).
"""
function animate_sensitivity(output_path)
    println("Rendering Sensitivity Animation to $output_path...")

    # Parameters for Chaos
    μ_val = -0.225
    ε_val = 0.25
    A_val = 0.3
    ω_val = 1.0

    p_vals = [μ_val, ε_val, A_val, ω_val]

    # Two close initial conditions
    u1_init = [0.1, 0.1]
    u2_init = [0.1 + 1e-6, 0.1]

    # Create CoupledODEs
    ds = get_ds(u1_init, p_vals)

    # We need to step manually to animate time evolution
    # Let's pre-compute the trajectories for performance
    T_total = 100.0
    dt = 0.05
    t_steps = 0:dt:T_total

    # Trajectory 1
    reinit!(ds, u1_init)
    tr1, t1 = trajectory(ds, T_total; Δt=dt, Ttr=0)

    # Trajectory 2
    reinit!(ds, u2_init)
    tr2, t2 = trajectory(ds, T_total; Δt=dt, Ttr=0)

    # Setup Figure
    fig = Figure(size=(1000, 500))
    ax1 = Axis(fig[1, 1], title="Phase Space", xlabel="x", ylabel="y")
    ax2 = Axis(fig[1, 2], title="Log Distance", xlabel="Time", ylabel="log10(|Δu|)")

    xlims!(ax1, -2.5, 2.5)
    ylims!(ax1, -2.0, 2.0)
    xlims!(ax2, 0, T_total)
    ylims!(ax2, -7, 1) # log10(1e-6) = -6, up to log10(approx 2) ~ 0.3

    # Observables for animation
    time_idx = Observable(1)

    # Tail length
    tail = 200

    # Plot tails
    points1 = @lift(Point2f.(tr1[max(1, $time_idx - tail):$time_idx, 1], tr1[max(1, $time_idx - tail):$time_idx, 2]))
    points2 = @lift(Point2f.(tr2[max(1, $time_idx - tail):$time_idx, 1], tr2[max(1, $time_idx - tail):$time_idx, 2]))

    lines!(ax1, points1, color=:cyan, linewidth=1.5, label="u1")
    lines!(ax1, points2, color=:magenta, linewidth=1.5, label="u2")

    # Plot heads
    head1 = @lift(Point2f(tr1[$time_idx, 1], tr1[$time_idx, 2]))
    head2 = @lift(Point2f(tr2[$time_idx, 1], tr2[$time_idx, 2]))

    scatter!(ax1, head1, color=:cyan, markersize=10)
    scatter!(ax1, head2, color=:magenta, markersize=10)

    # Plot Nullclines (Dynamic)
    nullcline_x = range(-2.5, 2.5, length=200)
    nc_points = @lift begin
        t_current = t_steps[$time_idx]
        pts = Point2f[]
        for x_val in nullcline_x
            denom = ε_val + μ_val * x_val^2
            if abs(denom) > 1e-4
                num = x_val - x_val^3 + A_val * cos(ω_val * t_current)
                y_val = num / denom
                if abs(y_val) < 100.0
                    push!(pts, Point2f(x_val, y_val))
                else
                    push!(pts, Point2f(x_val, NaN))
                end
            else
                push!(pts, Point2f(x_val, NaN))
            end
        end
        pts
    end

    # x-nullcline
    hlines!(ax1, [0.0], color=:orange, linestyle=:dash, linewidth=2, label="x-nullcline")
    # y-nullcline
    lines!(ax1, nc_points, color=:cyan, linestyle=:dash, linewidth=2, label="y-nullcline")

    # Add Legend
    axislegend(ax1, position=:rt)

    # Plot Distance
    # dist = sqrt((x1-x2)^2 + (y1-y2)^2)
    dists = [norm(tr1[i] - tr2[i]) for i in 1:min(length(tr1), length(tr2))]
    log_dists = log10.(dists .+ 1e-16) # Avoid log(0)

    # Plot full distance curve up to current time
    dist_points = @lift(Point2f.(t_steps[1:$time_idx], log_dists[1:$time_idx]))
    lines!(ax2, dist_points, color=:white, linewidth=1)

    # Record
    record(fig, output_path, 1:length(t_steps); framerate=30) do i
        time_idx[] = i
    end
    println("Saved $output_path")
end

"""
    plot_heatmap(A_vals, ω_vals, matrix, title_str, output_path)

Plots a heatmap of the chaos metric.
"""
function plot_heatmap(A_vals, ω_vals, matrix, title_str, output_path)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title=title_str, xlabel="Forcing Frequency (ω)", ylabel="Forcing Amplitude (A)")

    heatmap!(ax, ω_vals, A_vals, matrix, colormap=:inferno)
    Colorbar(fig[1, 2], limits=extrema(matrix), colormap=:inferno)

    save(output_path, fig)
    println("Saved $output_path")
end

end
