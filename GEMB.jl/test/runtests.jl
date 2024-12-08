using CondaPkg, PythonCall

using Kerchunk

# Raster creation and analysis packages
using Rasters, ArchGDAL # to save a raster
using YAXArrays # to open a raster

using Test

@testset "Kerchunk.jl" begin
    @static if !(Sys.iswindows())
        include("python_local_kerchunk.jl")
    end
    # include("its_live_catalog.jl")
    include("corrections.jl")
    include("real_zarr.jl")
    end
end
