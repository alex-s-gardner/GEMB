#=
# Solar Dynamics Observatory
=#

# First, we download the Kerchunk catalog:
using Downloads
Downloads.download("https://esip-qhub-public.s3-us-west-2.amazonaws.com/noaa/nwm/nwm_reanalysis.json.zst", "nwm_reanalysis.json.zst")

# This is a compressed file.  We can decompress it in memory using the [`TranscodingStreams.jl`](https://github.com/JuliaIO/TranscodingStreams.jl) API.
using CodecZstd
write("nwm_reanalysis.json", transcode(ZstdDecompressor, read("nwm_reanalysis.json.zst")))

# We can also open the decompressed catalog directly in Zarr.jl, or any package that sits on top of it.  Let's open it in plain Zarr first:
using Zarr
using ZarrDatasets, Rasters, YAXArrays
# We can open in Zarr directly, which gives us a ZarrGroup:
zg = Zarr.zopen("reference://nwm_reanalysis.json")
# or in ZarrDatasets.jl, which accounts for CF conventions:
zd = ZarrDataset("reference://nwm_reanalysis.json")
# Note the more descriptive display here.

# or in YAXArrays.jl, 
# ```julia
# ya = YAXArrays.open_dataset(...)
# ```
# And perform random operations
# Or visualize a subset (show this!!!!)