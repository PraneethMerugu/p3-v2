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
    dy = x - x^3 - ε * y - μ * x^2 * y + A * cos(ω * t)

    return SVector(dx, dy)
end
end
