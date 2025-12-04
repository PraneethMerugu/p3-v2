module Chaos
using DynamicalSystems
using ComplexityMeasures
using StateSpaceSets
using ..Model # Import the rule from above

using DifferentialEquations

export get_ds, compute_lz_complexity, compute_markov_matrix, compute_chaos_grid

"""
    get_ds(u0, p)

Returns a CoupledODEs (ContinuousDynamicalSystem) initialized with the Duffing rule.
Includes a safety callback to terminate if the system diverges (|u| > 100.0).
"""
function get_ds(u0=[0.1, 0.1], p=p_default)
    # 2. Safety Callback to prevent Finite-Time Blowup
    # Terminate if |x| or |y| exceeds a threshold (e.g., 100.0)
    is_unstable(u, t, integrator) = abs(u[1]) > 100.0 || abs(u[2]) > 100.0
    terminate_cb = DiscreteCallback(is_unstable, terminate!)

    # CoupledODEs is the standard struct for continuous systems in DynamicalSystems.jl v3+
    # We pass the callback via diffeq kwargs.
    # We also keep isoutofdomain as a fallback for functions like trajectory() that might ignore terminate!
    diffeq = (
        callback=terminate_cb,
        isoutofdomain=(u, p, t) -> abs(u[1]) > 1e4 || abs(u[2]) > 1e4
    )
    return CoupledODEs(duffing_rule, u0, p; diffeq)
end

"""
    compute_lz_complexity(ds, T; Ttr=1000.0)

Computes the Lempel-Ziv complexity of the system `ds` over time `T`.
Uses a simple binary symbolization based on the x-coordinate (x > 0).
"""
function compute_lz_complexity(ds::CoupledODEs, T::Real; Ttr::Real=1000.0)
    # 1. Generate trajectory
    # We need a discrete time series for LZ complexity.
    # Sampling period: approximate roughly 100 points per forcing period (2π/ω)
    # Default ω=1.0 => T_force ≈ 6.28. dt=0.05 is reasonable.
    dt = 0.05
    tr = nothing
    try
        tr, t = trajectory(ds, T; Δt=dt, Ttr=Ttr)
    catch e
        @warn "Trajectory generation failed (likely unstable): $e"
        return 0.0
    end

    # 2. Symbolize (Binary: x > 0)
    # We take the first dimension (x)
    x_series = tr[:, 1]
    # Binary symbolization: 1 if x > 0, 0 otherwise
    s = [x > 0 ? 1 : 0 for x in x_series]

    # 3. Compute LZ Complexity
    # Normalize by sequence length? The standard lempel_ziv_complexity function 
    # in ComplexityMeasures usually returns the raw count or normalized.
    # We'll use the standard function.
    c = complexity(LempelZiv76(), s)

    # Normalize (optional, but good for comparison)
    # C_norm = C / (N / log2(N))
    N = length(s)
    if N < 2
        return 0.0
    end
    norm_factor = N / log2(N)
    return c / norm_factor
end

"""
    compute_markov_matrix(ds, T; Ttr=1000.0, bins=20)

Computes a Markov transition matrix for the system.
Partitions the phase space into a grid and counts transitions.
"""
function compute_markov_matrix(ds::CoupledODEs, T::Real; Ttr::Real=1000.0, bins::Int=20)
    # 1. Generate trajectory
    dt = 0.1
    tr = nothing
    try
        tr, t = trajectory(ds, T; Δt=dt, Ttr=Ttr)
    catch e
        return zeros(Int, 0, 0)
    end

    # Check for NaNs or Infs
    if any(!isfinite, Matrix(tr))
        @warn "Trajectory contains non-finite values. Skipping Markov computation."
        return zeros(Int, 0, 0)
    end

    # 2. Partition Phase Space
    # We'll use a simple grid partition on the observed range
    # StateSpaceSets.RectangularBinning is efficient

    # Define binning
    # We want 'bins' subdivisions along each axis
    binning = ValueBinning(RectangularBinning(bins))

    # 3. Symbolize trajectory into visited bins
    # Use codify from ComplexityMeasures
    symbol_sequence = codify(binning, tr)

    # 4. Compute Transition Matrix
    # Identify unique states
    unique_states = unique(symbol_sequence)
    n_states = length(unique_states)
    state_map = Dict(s => i for (i, s) in enumerate(unique_states))

    M = zeros(Int, n_states, n_states)

    for i in 1:(length(symbol_sequence)-1)
        current_s = symbol_sequence[i]
        next_s = symbol_sequence[i+1]

        row = state_map[current_s]
        col = state_map[next_s]
        M[row, col] += 1
    end

    # Normalize rows to get probabilities
    P = zeros(Float64, n_states, n_states)
    for i in 1:n_states
        row_sum = sum(M[i, :])
        if row_sum > 0
            P[i, :] = M[i, :] ./ row_sum
        end
    end

    return P
end

"""
    compute_chaos_grid(duffing_sys, A_range, ω_range, fixed_p)

Performs a grid search over A and ω to map out chaos (using LZ complexity).
Note: duffing_sys is ignored here as we use the internal rule, but kept for signature compatibility if needed.
We'll just use the ranges.
"""
function compute_chaos_grid(duffing_sys, A_range, ω_range, fixed_p)
    # Initialize output matrix
    lz_matrix = zeros(length(A_range), length(ω_range))

    # Pre-allocate system to avoid overhead
    # Initial parameters (will be overwritten)
    p_init = [fixed_p[:μ], fixed_p[:ε], A_range[1], ω_range[1]]
    ds = get_ds([0.1, 0.1], p_init)

    # Indices in p vector
    # p = [μ, ε, A, ω]
    idx_A = 3
    idx_ω = 4

    # Loop
    # Parallelization could be added here with Threads.@threads
    # For now, simple loop
    for (j, ω_val) in enumerate(ω_range)
        for (i, A_val) in enumerate(A_range)
            # Update parameters
            set_parameter!(ds, idx_A, A_val)
            set_parameter!(ds, idx_ω, ω_val)

            # Re-initialize state to avoid getting stuck in far-away attractors?
            # Or continue from previous? Continuing is faster but might have hysteresis.
            # Let's re-init to a fixed point to be consistent.
            reinit!(ds, [0.1, 0.1])

            # Compute LZ
            # Short run for heatmap speed
            lz = compute_lz_complexity(ds, 500.0; Ttr=100.0)
            lz_matrix[i, j] = lz
        end
    end

    return lz_matrix
end

end
