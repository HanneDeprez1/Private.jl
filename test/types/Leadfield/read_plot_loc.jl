using Leadfield
using Winston
using TextPlots

positions, neighbours = readLOC("/Users/rluke/.julia/v0.3/Leadfield/data/FEM_Default50MRI.loc", verbose=true)

println()
println()
println("Back")
plot(vec(positions[1,:]), vec(positions[3,:]), cols=27, rows=15)


println()
println()
println("Side")
plot(vec(positions[2,:]), vec(positions[3,:]), cols=27, rows=15)


println()
println()
println("Top")
plot(vec(positions[1,:]), vec(positions[2,:]), cols=27, rows=15)
