#=
# NASA OISST data
=#

import AWSS3
AWS.global_aws_config(AWS.AWSConfig(; region = "us-east-1"))

catalog = download("https://ncsa.osn.xsede.org/Pangeo/pangeo-forge/pangeo-forge/aws-noaa-oisst-feedstock/aws-noaa-oisst-avhrr-only.zarr/reference.json")

using Kerchunk, Rasters, ZarrDatasets

rs = RasterStack("reference://$catalog"; source = Rasters.Zarrsource(), lazy = true)

# We can now manipulate the data using any Rasters or DimensionalData functions, or any Julia array manipulation functions