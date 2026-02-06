# `model_initialize_parameters` documentation
`model_initialize_parameters` initializes and validates the model configuration options, setting default values for physics modules, grid geometry, and output controls.

# Syntax
```matlab
ModelParam = model_initialize_parameters
ModelParam = model_initialize_parameters(option1=value1, ..., optionN=valueN)
```

# Description

`ModelParam = model_initialize_parameters` creates a `ModelParam` structure containing 38 fields of default model parameters. 

`ModelParam = model_initialize_parameters(option1=value1,...,optionN=valueN)` creates a `ModelParam` structure with specified parameters set to preferred values. Any unspecified parameters are set to default values. 

# Inputs 
A complete list of GEMB variables and descriptions can be found [here](GEMB_variables.md). The fields pertinent to `model_initialize_parameters` are described below:

```matlab
--- GENERAL & INITIALIZATION ---
.run_prefix                    : string       Unique run identifier (default: "default").
.spinup_cycles                 : integer      Number of spin-up cycles (default: 0).

--- DENSITY & DENSIFICATION ---
.densification_method          : string       Model: "HerronLangway", "Anthern", "Ligtenberg".
.densification_coeffs_M01      : string       Coeffs for Ligtenberg model (e.g., "Gre_RACMO_GS_SW0").
.new_snow_method               : string       Fresh snow density model (e.g., "350kgm2", "Fausto").
.density_ice                   : double       Density of glacier ice [kg m^-3].
.rain_temperature_threshold    : double       Temperature threshold [C] above which precipitation falls as rain.

--- LONGWAVE EMISSIVITY & ROUGHNESS ---
.emissivity_method             : string       Method: "uniform", "re_threshold", "re_w_threshold".
.emissivity                    : double       Base longwave emissivity (0-1).
.emissivity_grain_grain_radius_large   : double       Emissivity for large grain sizes (0-1).
.emissivity_grain_radius_threshold: double       Grain radius threshold [mm] for emissivity switching.
.surface_roughness_effective_ratio : double   Ratio of physical surface roughness to effective roughness.

--- THERMAL CONDUCTIVITY ---
.thermal_conductivity_method   : string       Model: "Sturm" or "Calonne".

--- MELT & LIQUID WATER ---
.water_irreducible_saturation  : double       Irreducible water content saturation fraction (0-0.2).

--- ALBEDO & RADIATION ---
.albedo_method                 : string       Scheme: "GardnerSharp", "GreuellKonzelmann", etc.
.albedo_density_threshold      : double       Density threshold [kg m^-3] below which albedo_method is applied (Default Inf).
.shortwave_absorption_method   : double       0 (surface only) or 1 (subsurface penetration).
.albedo_snow                   : double       Albedo for fresh snow (0.5-0.95).
.albedo_ice                    : double       Albedo for bare ice (0.2-0.6).
.albedo_fixed                  : double       Fixed albedo used if albedo_method="None" or density > threshold.
.shortwave_downward_diffuse    : double       Downward diffusive shortwave flux [W m^-2].
.solar_zenith_angle            : double       Solar Zenith Angle [degrees].
.cloud_optical_thickness       : double       Cloud Optical Thickness.
.black_carbon_snow             : double       BC concentration in snow [ppm].
.black_carbon_ice              : double       BC concentration in ice [ppm].
.cloud_fraction                : double       Cloud fraction (0-1) for GreuellKonzelmann method.
.albedo_wet_snow_t0            : double       Time scale [d] for wet snow (Bougamont2005).
.albedo_dry_snow_t0            : double       Time scale [d] for dry snow (Bougamont2005).
.albedo_K                      : double       Temperature coef. time scale [d] (Bougamont2005).

--- OUTPUT CONTROLS ---
.output_frequency              : string       Output resolution: "daily", "monthly", or "all".
.output_padding                : integer      Extra vertical levels for grid resizing (default: 1000).

--- GRID GEOMETRY ---
.column_ztop                   : double       Depth of constant grid spacing at the surface [m].
.column_dztop                  : double       Initial surface grid spacing [m].
.column_dzmin                  : double       Minimum allowable grid spacing [m].
.column_dzmax                  : double       Maximum allowable grid spacing [m].
.column_zmax                   : double       Maximum total column depth [m].
.column_zmin                   : double       Minimum total column depth [m].
.column_zy                     : double       Grid stretching factor for lower layers.
```

# Column Geometry 

![](https://gmd.copernicus.org/articles/16/2277/2023/gmd-16-2277-2023-f01-web.png)

Grid geometry is set by the `column_*` inputs, following the variable names in Figure 1 of [*Gardner et al.,* 2023](https://doi.org/10.5194/gmd-16-2277-2023).


# Examples 
## Example 1: Default values
Create a `ModelParam` structure to containing default values: 

```matlab
>> ModelParam = model_initialize_parameters
ModelParam = 
  struct with fields:

                           run_prefix: "default"
                        spinup_cycles: 0
                 densification_method: "Anthern"
             densification_coeffs_M01: "Gre_RACMO_GS_SW0"
                      new_snow_method: "350kgm2"
                          density_ice: 910.00
           rain_temperature_threshold: 273.15
                    emissivity_method: "uniform"
    surface_roughness_effective_ratio: 0.10
                           emissivity: 0.97
        emissivity_grain_radius_large: 0.97
    emissivity_grain_radius_threshold: 10.00
          thermal_conductivity_method: "Sturm"
         water_irreducible_saturation: 0.07
                        albedo_method: "GardnerSharp"
             albedo_density_threshold: Inf
          shortwave_absorption_method: 0
                          albedo_snow: 0.85
                           albedo_ice: 0.48
                         albedo_fixed: 0.85
           shortwave_downward_diffuse: 0
                   solar_zenith_angle: 0
              cloud_optical_thickness: 0
                    black_carbon_snow: 0
                     black_carbon_ice: 0
                       cloud_fraction: 0.10
                   albedo_wet_snow_t0: 15.00
                   albedo_dry_snow_t0: 30.00
                             albedo_K: 7.00
                     output_frequency: "all"
                       output_padding: 1000.00
                          column_ztop: 10.00
                         column_dztop: 0.05
                         column_dzmin: 0.03
                         column_dzmax: 0.08
                          column_zmax: 250.00
                          column_zmin: 130.00
                            column_zy: 1.10
```

## Example 2: Override defaults
If you wish to override any default model parameters, you can specify them as keyword arguments. For example, here's how to specify three spinup cycles instead of the default zero: 

```matlab
ModelParam = model_initialize_parameters(spinup_cycles=3);
```
To override multiple default values, list them as function inputs like this: 

```matlab
ModelParam = model_initialize_parameters(spinup_cycles=3,...
    ice_density=920,...
    densification_method="Ligtenberg");
```
## Example 3: Update model parameters
If you love your `ModelParam` structure very much, but you just want to tweak one little thing about it, you can update any field value manually. Just make sure spell and capitalize the field name correctly. Here's how to initialize a `ModelParam` structure, then change the `column_ztop` value from the default 10 to 15: 

```matlab
% Initialize a model parameter structure: 
ModelParam = model_initialize_parameters;

% Override the default column_ztop value: 
ModelParam.column_ztop = 15;
```

# Accessing help
Not sure what your options are? Try entering "help" as the selection: 

```matlab
model_initialize_parameters(densification_method="options")
Error using model_initialize_parameters (line 131)
model_initialize_parameters(densification_method="options")

Invalid value for 'densification_method' argument. Value must be a member of this set:
'HerronLangway'
'Anthern'
'Ligtenberg' 
```
Above, you see there's nothing special about the word "help" per se, but it does produce an error message containing a helpful list of `densification_method` options. 


# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
