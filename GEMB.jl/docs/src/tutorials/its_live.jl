#=
# ITS_LIVE data

This package's development was motivated by the [ITS_LIVE](https://its-live.jpl.nasa.gov/) project.
=#

using Rasters       # Raster data analysis in Julia
using ZarrDatasets  # Zarr support for Rasters
using Kerchunk      # Kerchunk support for Zarr
using Statistics    # Basic statistics
# We can load the catalog from the catalog file that we ship in the repo:
catalog_path = joinpath(dirname(dirname(pathof(Kerchunk))), "test", "data", "its_live_catalog.json")
rs = RasterStack("reference://$(catalog_path)"; source = Rasters.Zarrsource())
# We've now loaded the dataset lazily in a `RasterStack`, which is essentially a stack of multiple variables.
# Now, we can apply arbitrary Rasters.jl functions to the stack, or plot it, and treat it as a general Julia array!
# 
# Let's plot first:
using CairoMakie
heatmap(rs.v)
# We can also aggregate the data to a lower resolution, which downloads the entire dataset.
# Here, we aggregate by a factor of 10 in both dimensions, so a 10x10 window is aggregated to a single pixel.
vs2 = Rasters.aggregate(rs, mean, 10) # now everything is loaded in disk
# and plot this aggregated data:
arrows(dims(vs2, X) |> collect, dims(vs2, Y) |> collect, vs2.vx .* 20, vs2.vy .* 20; arrowsize = 5)