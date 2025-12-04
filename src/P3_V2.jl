module P3_V2

# Core Model
include("model.jl")
using .Model
export duffing_rule, p_default

# Modules
include("stability.jl")
using .Stability
export get_jacobian_at_point, analyze_equilibrium

include("chaos.jl")
using .Chaos
export get_ds

# include("utils.jl")
# using .Utils
# export set_project_theme!

# include("visualizations.jl")
# using .Visualizations
# export animate_collision, animate_sensitivity, plot_heatmap

end
