```@meta
CurrentModule = Kerchunk
```

# Kerchunk

Kerchunk.jl is a Julia package that enables loading [Kerchunk reference catalogs](https://fsspec.github.io/kerchunk/) as [Zarr.jl](https://github.com/JuliaIO/Zarr.jl) arrays.


## Quick start

Kerchunk.jl is simply a storage backend to [`Zarr.jl`](https://github.com/JuliaIO/Zarr.jl).  Zarr does integrate with the more fully featured packages [`Rasters.jl`](https://github.com/rafaqz/Rasters.jl) and [`YAXArrays.jl`](https://github.com/JuliaDataCubes/YAXArrays.jl), which are the packages you will want to use to interact with Kerchunk data.

```julia
using Kerchunk, Zarr

za = Zarr.zopen("reference://path/to/kerchunk/catalog.json")
# and treat it like any other Zarr array!
# You can even wrap it in YAXArrays.jl to get DimensionalData.jl accessors:
using YAXArrays
YAXArrays.open_dataset(za)
# or open it as a Rasters.RasterStack:
using Rasters
Rasters.RasterStack(
    "reference://catalog.json", 
    source = Rasters.Zarrsource(),
    lazy = true, # need to include this
) # source must be explicit
```

It's most useful to open Kerchunk datasets as either RasterStacks or YAXArrays datasets, since both of those packages have great dimensionality support.

## Background

[`kerchunk`](https://fsspec.github.io/kerchunk/) is a Python package that generates the reference catalogs.

## Limitations
- No support for `gen` references with templates.
- No support for complex Jinja2 templates in `refs`.  (Although Kerchunk hardly supports this either...)

## Acknowledgements

This effort was funded by the NASA MEaSUREs program in contribution to the Inter-mission Time Series of Land Ice Velocity and Elevation (ITS_LIVE) project (https://its-live.jpl.nasa.gov/).

## Alternatives and related packages

- You can always use Python's `xarray` directly via PythonCall.jl
- [FSSpec.jl](https://github.com/asinghvi17/FSSpec.jl) is an alternative storage backend for Zarr.jl that wraps the same [`fsspec`](https://github.com/fsspec/filesystem_spec) that `xarray` uses under the hood.

This package is of course built on top of [Zarr.jl](https://github.com/JuliaIO/Zarr.jl), which is a pure-Julia Zarr array library.
[YAXArrays.jl](https://github.com/JuliaDataCubes/YAXArrays.jl) is a Julia package that can wrap Zarr arrays in a DimensionalData-compatible interface.
