# `model_initialize_profile` documentation
`model_initialize_profile ` initializes a GEMB column based on specified model and climate forcing parameters.

See also [`gemb_profile`](gemb_profile_documentation.md) and [`gemb_spinup`](gemb_spinup_documentation.md) 

# Syntax 
```matlab
Profile = model_initialize_profile(ModelParam, ClimateForcing)
```

# Description
`Profile = model_initialize_profile(ModelParam, ClimateForcing)` uses inputs `ModelParam` from [`model_initialize_parameters`](model_initialize_parameters_documentation.md) and input timetable `ClimateForcing` (which must contain at least `temperature_air_mean`). The output `Profile` is a table containing the following variables:  

|Variable   |  Units | Description          |      Initial Value|
|---|---|---|---|
|`z_center`    |   m  | Center height of each grid cell |                  array generated from model parameters |
|`dz`    |   m  |grid cell thickness |                  array generated from model parameters |
|`temperature`     | K |  temperature |  mean annual surface temperature |
|`density`     | kg m<sup>-3</sup> | ice density   |             `ModelParam.density_ice` |
|`water`     | kg m<sup>-2</sup> | old-snow water content  |   0|
|`grain_radius`    | mm  |old-snow grain size   |     2.5|
|`grain_dendricity` | fraction | old-snow grain dendricity | 1|
|`grain_sphericity` | fraction | old-snow grain sphericity|  0.5|
|`albedo` | fraction | surface albedo         |    `ModelParam.albedo_snow`|
|`albedo_diffuse` | fraction | diffuse surface albedo  |   `ModelParam.albedo_snow` |

# Example
Initialize a GEMB column based on default model parameters and a mean surface temperature of -20 C. 

```matlab
% Initialize parameters: 
ModelParam = model_initialize_parameters;

% Mean air temperature is the only ClimateForcing parameter required to create a Profile:
ClimateForcing.Properties.CustomProperties.temperature_air_mean = 253.15; % -20 C

% Initialize Column: 
Profile = model_initialize_profile(ModelParam, ClimateForcing);
```

Use MATLAB's built-in `head` function to view the top rows of the `Profile` table:

```matlab
>> head(Profile)
     z_center     dz     temperature    density    water    grain_radius    grain_dendricity    grain_sphericity    albedo    albedo_diffuse
     ________    ____    ___________    _______    _____    ____________    ________________    ________________    ______    ______________
      -0.025     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.075     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.125     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.175     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.225     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.275     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.325     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85     
      -0.375     0.05      253.15         910        0          2.5                1                  0.5            0.85          0.85               
```

Above, we see that all ten variables are 264x1, representing a column of initial values. Here is how the initial grid-cell thickness `dz` and `temperature` values vary as a function of depth:

```matlab
figure
subplot(1,2,1)
plot(Profile.dz,Profile.z_center,'o-')
xlabel 'Vertical spacing (m)'
ylabel 'Grid cell center height (m)'

subplot(1,2,2)
plot(Profile.temperature,Profile.z_center,'o-')
xlabel 'Initial temperature (K)'

set(gcf,'Renderer','painters')
exportgraphics(gcf,'model_initialize_profile_documentation_01.jpg','resolution',500)
```

![](https://chadagreene.com/GEMB_figures/model_initialize_profile_documentation_01.jpg)


# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
