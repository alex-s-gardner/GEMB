# module GEMBDirect

using Dates, Statistics
using ProgressMeter

include("accumulation.jl")
include("albedo.jl")
include("densification.jl")
include("gardner_alb.jl")
include("grain_growth.jl")
include("grid_initialize.jl")
include("managelayers.jl")
include("melt.jl")
include("shortwave.jl")
include("thermo.jl")

include("GEMB.jl")



# end