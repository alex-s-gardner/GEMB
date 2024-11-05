#=
# Creating Kerchunk catalogs

Kerchunk.jl is only a Kerchunk reader, meaning that if you want to create Kerchunk catalogs, you need to use the Python `kerchunk` package.

The easiest way to do this in Julia is to use the [CondaPkg.jl](https://github.com/JuliaPy/CondaPkg.jl) package to install the `kerchunk` package into a Conda environment and then use PythonCall.jl to call the `kerchunk` package.  This ensures reproducibility, since you can pin the versions in the generated CondaPkg.toml as well, and package management via CondaPkg has a very similar interface to `Pkg.jl`.
=#


#=

## Setting up the Conda environment

```julia
using CondaPkg
CondaPkg.add("python")
CondaPkg.add("kerchunk")
```


=#

#=

## Creating a Kerchunk catalog

=#

using CondaPkg, PythonCall

# There are two approaches to this - either call Python explicitly via the command line, or call Python via PythonCall.jl.
# Calling Python via Julia is a nicer interface, but you will very quickly run into binary incompatibility issues.

CondaPkg.withenv() do
run(```
$(CondaPkg.which("python")) -e "
import kerchunk
import kerchunk.hdf5 as hdf

# do something
"
```)
end

#=

## Using PythonCall.jl

```julia

```

=#

# Let's load this using Kerchunk.jl, and see what we get!

using Zarr, Kerchunk
z = Zarr.zopen("reference://catalog.json")

using Rasters, ZarrDatasets
r = Raster(z)