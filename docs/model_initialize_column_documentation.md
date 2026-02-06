# `model_initialize_column` documentation
`model_initialize_column ` initializes a GEMB column based on specifiedmodel and climate forcing parameters.

# Syntax 
```matlab
[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
```

# Description
`[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing)` uses inputs `ModelParam` from `model_initialize_parameters` and input structure `ClimateForcing` containing the field `temperature_air_mean`. Outputs are as follows: 

|Variable   |  Units | Description          |      Initial Value|
|---|---|---|---|
|`temperature`     | K |  temperature |  mean annual surface temperature |
|`dz`    |   m  |thickness |                  array generated from model parameters |
|`density`     | kg m^-3 | ice density   |             `ModelParam.density_ice` |
|`water`     | kg m^-2 | old-snow water content  |   0|
|`grain_radius`    | mm  |old-snow grain size   |     2.5|
|`grain_dendricity` | fraction | old-snow grain dendricity | 0|
|`grain_sphericity` | fraction | old-snow grain sphericity|  0|
|`albedo` | fraction | surface albedo         |    `ModelParam.albedo_snow`|
|`albedo_diffuse` | fraction | diffuse surface albedo  |   `ModelParam.albedo_snow` |

# Example
Initialize a GEMB column based on default model parameters and a mean surface temperature of -20 C. 

```matlab
% Initialize parameters: 
ModelParam = model_initialize_parameters;
ClimateForcing.temperature_air_mean = 253.15; % -20 C

% Initialize Column: 
[temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
```

Use MATLAB's built-in `whos` function to sizes of the arrays that were created by `model_initialize_column`:

```matlab
>> whos
  Name                    Size            Bytes  Class     Attributes

  ClimateForcing          1x1               176  struct              
  ModelParam              1x1              7972  struct              
  albedo                264x1              2112  double              
  albedo_diffuse        264x1              2112  double              
  density               264x1              2112  double              
  dz                    264x1              2112  double              
  grain_dendricity      264x1              2112  double              
  grain_radius          264x1              2112  double              
  grain_sphericity      264x1              2112  double              
  temperature           264x1              2112  double              
  water                 264x1              2112  double                         
```

Above, we see that all of the variables creates by `model_initialize_column` are 264x1, representing a column of initial values. Use [`dz2z`](dz2z_documentation.md) to convert the column of `dz` values to grid cell centers and show the grid spacing alongside initial temperature:

```matlab
z = dz2z(dz); 

figure
subplot(1,2,1)
plot(dz,z,'o-')
xlabel 'Vertical spacing (m)'
ylabel 'Grid cell center height (m)'

subplot(1,2,2)
plot(T,z,'o-')
xlabel 'Initial temperature (K)'

set(gcf,'Renderer','painters')
exportgraphics(gcf,'model_initialize_column_documentation_01.jpg','resolution',500)
```

![](https://chadagreene.com/GEMB_figures/model_initialize_column_documentation_01.jpg)


# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
