# ERA5 time series data 

Some of the GEMB documentation uses ERA5 Land hourly time-series data from a location near Greenland's [Summit Station](https://en.wikipedia.org/wiki/Summit_Camp). You can grab the exact same data files that are used in the GEMB documtation [here](https://chadagreene.com/GEMB_data/), but if you'd like to explore other locations, here's where and how we got the example data. 

The first step is to go to the ERA time series data page here: [https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries?tab=download](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries?tab=download)

![](https://chadagreene.com/GEMB_figures/ERA5_download_screenshot.png)

After you've logged in, click on the **Download** tab and select which variables to download. All of the following variables are necessary to run GEMB except skin temperature, which is only used for comparison to GEMB output. 

* 2m dewpoint temperature 
* 2m temperature
* Surface pressure
* Total precipitation (de-accumulated)
* Surface solar radiation downwards (de-accumulated)
* Surface thermal radiation downwards (de-accumulated)
* 10m u-component of wind
* 10m v-component of wind
* Skin temperature (optional)

Summit Station is located at (72.579583°N, 38.459186°W), which we round to the nearest tenth of a degree to get the closest ERA5 grid cell. For our example data, we selected hourly data from the start of the year 2000 through the end of 2025. Select NetCDF as the output format and request to download the data. The request may take a minute or two, but should be ready for download in less time than it takes to refill a coffee. 

After downloading, unzip the folder and you will find [five NetCDF files](https://chadagreene.com/GEMB_data/)with names like `reanalysis-era5-land-timeseries-sfc-2m-temperaturelaos9yiu.nc`. The last eight characters are unique identfiers automatically assigned to each file. Despite requesting all of the files at the same time, the unique identifiers are different for each of the five files, which can make it difficult to identify them as belonging together. To make things easier, we manually changed the final characters of all five filenames to `_summit`. 

In MATLAB, you can explore the contents of each NetCDF using the `ncdisp` function like this: 

```matlab
>> ncdisp('reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc')
Source:
           /Users/cgreene/Documents/GEMB_data/reanalysis-era5-land-timeseries-sfc-2m-temperaturelaos9yiu.nc
Format:
           netcdf4
Global Attributes:
           Conventions            = 'CF-1.7'
           GRIB_centre            = 'ecmf'
           GRIB_centreDescription = 'European Centre for Medium-Range Weather Forecasts'
           GRIB_edition           = 1
           GRIB_subCentre         = 0
           history                = '2025-02-10T00:00 GRIB to CDM+CF via cfgrib-0.9.10.4/ecCodes-2.26.0 with {"source": "../../../../../admp-scratch/data/arco-job-dcc56846-742c-4e86-aa97-2a9ae4af5f6d-inputs/pre-reanalysis_era5_land-sfc-2m-temperature:0-t2m-20250205-000000.grib", "filter_by_keys": {"stepType": "instant"}, "encode_cf": ["parameter", "time", "geography", "vertical"]}'
           institution            = 'European Centre for Medium-Range Weather Forecasts'
Dimensions:
           valid_time = 227928
Variables:
    d2m       
           Size:       227928x1
           Dimensions: valid_time
           Datatype:   single
           Attributes:
                       _FillValue                              = NaN
                       GRIB_NV                                 = 0
                       GRIB_Nx                                 = 3600
                       GRIB_Ny                                 = 1801
                       GRIB_cfName                             = 'unknown'
                       GRIB_cfVarName                          = 'd2m'
                       GRIB_dataType                           = 'fc'
                       GRIB_gridDefinitionDescription          = 'Latitude/Longitude Grid'
                       GRIB_gridType                           = 'regular_ll'
                       GRIB_iDirectionIncrementInDegrees       = 0.1
                       GRIB_iScansNegatively                   = 0
                       GRIB_jDirectionIncrementInDegrees       = 0.1
                       GRIB_jPointsAreConsecutive              = 0
                       GRIB_jScansPositively                   = 0
                       GRIB_latitudeOfFirstGridPointInDegrees  = 90
                       GRIB_latitudeOfLastGridPointInDegrees   = -90
                       GRIB_longitudeOfFirstGridPointInDegrees = 0
                       GRIB_longitudeOfLastGridPointInDegrees  = 359.9
                       GRIB_missingValue                       = 3.402823466385289e+38
                       GRIB_name                               = '2 metre dewpoint temperature'
                       GRIB_numberOfPoints                     = 6483600
                       GRIB_paramId                            = 168
                       GRIB_shortName                          = '2d'
                       GRIB_stepType                           = 'instant'
                       GRIB_stepUnits                          = 1
                       GRIB_totalNumber                        = 0
                       GRIB_typeOfLevel                        = 'surface'
                       GRIB_units                              = 'K'
                       long_name                               = '2 metre dewpoint temperature'
                       standard_name                           = 'unknown'
                       units                                   = 'K'
                       coordinates                             = 'latitude longitude'
    t2m       
           Size:       227928x1
           Dimensions: valid_time
           Datatype:   single
           Attributes:
                       _FillValue                              = NaN
                       GRIB_NV                                 = 0
                       GRIB_Nx                                 = 3600
                       GRIB_Ny                                 = 1801
                       GRIB_cfName                             = 'unknown'
                       GRIB_cfVarName                          = 't2m'
                       GRIB_dataType                           = 'fc'
                       GRIB_gridDefinitionDescription          = 'Latitude/Longitude Grid'
                       GRIB_gridType                           = 'regular_ll'
                       GRIB_iDirectionIncrementInDegrees       = 0.1
                       GRIB_iScansNegatively                   = 0
                       GRIB_jDirectionIncrementInDegrees       = 0.1
                       GRIB_jPointsAreConsecutive              = 0
                       GRIB_jScansPositively                   = 0
                       GRIB_latitudeOfFirstGridPointInDegrees  = 90
                       GRIB_latitudeOfLastGridPointInDegrees   = -90
                       GRIB_longitudeOfFirstGridPointInDegrees = 0
                       GRIB_longitudeOfLastGridPointInDegrees  = 359.9
                       GRIB_missingValue                       = 3.402823466385289e+38
                       GRIB_name                               = '2 metre temperature'
                       GRIB_numberOfPoints                     = 6483600
                       GRIB_paramId                            = 167
                       GRIB_shortName                          = '2t'
                       GRIB_stepType                           = 'instant'
                       GRIB_stepUnits                          = 1
                       GRIB_totalNumber                        = 0
                       GRIB_typeOfLevel                        = 'surface'
                       GRIB_units                              = 'K'
                       long_name                               = '2 metre temperature'
                       standard_name                           = 'unknown'
                       units                                   = 'K'
                       coordinates                             = 'latitude longitude'
    latitude  
           Size:       1x1
           Dimensions: 
           Datatype:   double
           Attributes:
                       _FillValue       = NaN
                       long_name        = 'latitude'
                       standard_name    = 'latitude'
                       stored_direction = 'decreasing'
                       units            = 'degrees_north'
    longitude 
           Size:       1x1
           Dimensions: 
           Datatype:   double
           Attributes:
                       _FillValue    = NaN
                       long_name     = 'longitude'
                       standard_name = 'longitude'
                       units         = 'degrees_east'
    valid_time
           Size:       227928x1
           Dimensions: valid_time
           Datatype:   int64
           Attributes:
                       units    = 'hours since 1970-01-01'
                       calendar = 'proleptic_gregorian'
>> 
```

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 