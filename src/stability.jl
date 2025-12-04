module Stability
using LinearAlgebra
using ForwardDiff
using ..Model

export get_jacobian_at_point, analyze_equilibrium

# 1. Function to get numerical Jacobian at a specific point
function get_jacobian_at_point(u, p, t=0.0)
    # u: [x, y]
    # p: [μ, ε, A, ω]

    # Define a wrapper for ForwardDiff
    f = (u_val) -> duffing_rule(u_val, p, t)

    # Calculate Jacobian using ForwardDiff
    J = ForwardDiff.jacobian(f, u)

    return J
end

# 2. Analyze stability of a point
function analyze_equilibrium(u, p, t=0.0)
    J = get_jacobian_at_point(u, p, t)
    evals = eigen(J).values

    # Classification
    real_parts = real.(evals)
    imag_parts = imag.(evals)

    classification = "Unknown"
    if all(r -> r < 0, real_parts)
        classification = "Stable Sink"
        if any(i -> i != 0, imag_parts)
            classification = "Stable Spiral"
        end
    elseif all(r -> r > 0, real_parts)
        classification = "Unstable Source"
        if any(i -> i != 0, imag_parts)
            classification = "Unstable Spiral"
        end
    elseif any(r -> r > 0, real_parts) && any(r -> r < 0, real_parts)
        classification = "Saddle"
    elseif any(r -> r == 0, real_parts)
        classification = "Non-hyperbolic (Bifurcation Candidate)"
    end

    return (eigenvalues=evals, classification=classification, jacobian=J)
end
end
