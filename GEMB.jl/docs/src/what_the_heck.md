# What is Kerchunk?

Kerchunk is a powerful tool designed to optimize access to large scientific datasets, particularly those stored in cloud-based object stores. It addresses the challenges of working with numerous small files or large, chunked files by creating a unified, efficient interface for data access.

At its core, Kerchunk works by creating a "fake" file system that maps to a Zarr store.  The file system describes a mapping from Zarr chunks to byte ranges of the source files.  

This approach allows Kerchunk to effectively wrap one or many data files into a single Zarr array, providing a consolidated view of the data. By doing so, it enables faster data access, reduces the number of API calls needed to retrieve information (by essentially front-loading the process), and greatly simplifies the process of working with multi-file datasets.


## Available data sources

The unit of Kerchunking is the _catalog_.  Each catalog is either a single JSON file or a directory of Parquet files.  The catalog is essentially a dictionary of file paths mapped to byte ranges. 

Catalogs are "sidecar" files, and may not always be present with the original data.  Generally, at least for now, if there's no obvious Kerchunk file you would have to generate one yourself, using Python.  Julia support for constructing catalogs is currently nonexistent, but it's on the bucket list!

## Tips and tricks

### Where's my CRS?

That's an interesting question.  Over the short term, Julia doesn't have support for CF-style (climate-and-forecast conventions) CRS metadata.  Additionally, CRS from e.g NetCDF files are stored as empty variables, which Kerchunk removes.  

There are two places you might look for CRS information. 
- First, see if the global metadata contains a `crs_wkt` or `spatial_ref` field.  If so, you can use that.  Other potential keys to look for are `proj4string`, `proj4text`, or `spatial_epsg`.
- Second, you might find a `grid_mapping` metadata field in a layer / Zarr group, which will contain a link to the CRS.  If the value of that `grid_mapping` field is also a key in the global metadata, then that will contain the CRS.

If you're using `Rasters.jl` to load the data, you can set the CRS on a `Raster` or `RasterStack` like this:

```julia
ras = Rasters.setcrs(ras, new_crs)
```

and if you have a WKT string, for example, you can construct `new_crs` like this:

```julia
new_crs = Rasters.ESRIWellKnownText(wkt_string)
# or
new_crs = Rasters.EPSG(epsg_code)
# or
new_crs = Rasters.ProjString(proj4_string)
```

### S3 redirect errors

Many S3 buckets are restricted to only allow access from certain regions.  If you get an error like this:
```
nested task error: AWS.AWSExceptions.AWSException: PermanentRedirect -- The bucket you are attempting to access must be addressed using the specified endpoint. Please send all future requests to this endpoint.

HTTP.Exceptions.StatusError(301, "GET", "/its-live-data/velocity_image_pair/landsatOLI/v02/N70W040/LC08_L1GT_004010_20140206_20200912_02_T2_X_LC08_L1GT_004010_20140529_20200911_02_T2_G0120V02_P008.nc", HTTP.Messages.Response:
"""
HTTP/1.1 301 Moved Permanently
...
```

then you can set your AWS config to say you're coming from a different region, like this:
```
import AWS
AWS.global_aws_config(AWS.AWSConfig(; region="us-west-2"))
```

### Version mismatches between Python and Julia

Python and Julia load different versions of libraries, which can cause incompatibilities.  For example, both NCDatasets.jl and Python's netcdf4 library depend on libhdf5, but the versions they try to load are incompatible.

So, we recommend invoking Python via the shell and `run`, rather than using PyCall directly to Kerchunk datasets.