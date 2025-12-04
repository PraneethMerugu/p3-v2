using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using P3_V2
using GLMakie
using DynamicalSystems

println("Verifying Vector Field...")

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], title="Vector Field Verification")

# Test parameters
μ = -0.225
ε = Observable(0.25)
A = Observable(0.3)
ω = 1.0

# Call the function
plot_vector_field!(ax, μ, ε, A, ω)

# Add a slider to test reactivity
sg = SliderGrid(fig[2, 1],
    (label="A", range=0:0.01:1.0, startvalue=0.3)
)
sl = sg.sliders[1]
connect!(A, sl.value)

# Save a frame
save(joinpath(@__DIR__, "../verify_vector_field.png"), fig)
println("Verification successful. Saved verify_vector_field.png")
