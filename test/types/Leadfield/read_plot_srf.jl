using Leadfield
using Winston

Xcords, Ycords, Zcords = readSRF("/Users/rluke/Documents/GITs/Leadfield.jl/data/Default50Brain.srf");

p = plotSRF(Xcords, Ycords, Zcords);
file(p, "Brain.eps", width=1200, height=600)


Xcords, Ycords, Zcords = readSRF("/Users/rluke/Documents/GITs/Leadfield.jl/data/Default50Skin.srf");

p = plotSRF(Xcords, Ycords, Zcords);
file(p, "Skin.eps", width=1200, height=600)


Xcords, Ycords, Zcords = readSRF("/Users/rluke/Documents/GITs/Leadfield.jl/data/FEM_Default50MRI.srf");

p = plotSRF(Xcords, Ycords, Zcords);
file(p, "UK.eps", width=1200, height=600)
