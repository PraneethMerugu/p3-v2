module Chaos
using DynamicalSystems
using ..Model # Import the rule from above

export get_ds

"""
get_ds(u0, p)

Returns a CoupledODEs (ContinuousDynamicalSystem) initialized with the Duffing rule.
"""
function get_ds(u0=[0.1, 0.1], p=p_default)
    # CoupledODEs is the standard struct for continuous systems in DynamicalSystems.jl v3+
    return CoupledODEs(duffing_rule, u0, p)
end
end
