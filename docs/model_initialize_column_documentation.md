# `model_initialize_column` documentation
`model_initialize_column ` initializes a GEMB column based on specifiedmodel and climate forcing parameters.

# Syntax 
```matlab
[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing)
```

# Description
`[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing)` uses inputs `ModelParam` from `model_initialize_parameters` and input structure `ClimateForcing` containing the field `T_air_mean`. Outputs are as follows: 

|Variable   |  Units | Description          |      Initial Value|
|---|---|---|---|
|`T`     | K |  temperature |  mean annual surface temperature |
|`dz`    |   m  |thickness |                  array generated from model parameters |
|`d`     | kg m^-3 | ice density   |             `ModelParam.density_ice` |
|`W`     | kg m^-2 | old-snow water content  |   0|
|`re`    | mm  |old-snow grain size   |     2.5|
|`gdn` | fraction | old-snow grain dendricity | 0|
|`gsp` | fraction | old-snow grain sphericity|  0|
|`a` | fraction | surface albedo         |    `ModelParam.albedo_snow`|
|`a_diffuse` | fraction | diffuse surface albedo  |   `ModelParam.albedo_snow` |


# Example
Initialize a GEMB column based on default model parameters and a mean surface temperature of -20 C. 

```matlab
% Initialize parameters: 
ModelParam = model_initialize_parameters();
ClimateForcing.T_air_mean = 253.15; % -20 C

% Initialize Column: 
[T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing);
```

Use MATLAB's built-in `whos` function to sizes of the arrays that were created by `model_initialize_column`:

```matlab
>> whos
Name                  Size            Bytes  Class     Attributes

ClimateForcing        1x1               176  struct              
ModelParam            1x1              7972  struct              
T                   264x1              2112  double              
W                   264x1              2112  double              
a                   264x1              2112  double              
a_diffuse           264x1              2112  double              
d                   264x1              2112  double              
dz                  264x1              2112  double              
gdn                 264x1              2112  double              
gsp                 264x1              2112  double              
re                  264x1              2112  double              
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
