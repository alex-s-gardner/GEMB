# `gemb_interp` documentation
`gemb_interp`  regularizes the MxN gemb output onto a consistent vertical grid. This function is necessary because the vertical spacing of GEMB output evolves with every timestep. 

# Syntax 
```matlab
A_regularized = gemb_interp(Z_center,A,Profile)
A_regularized = gemb_interp(Z_center,A,Profile,interp_method)
```

# Description

`A_regularized = gemb_interp(Z_center,A,Profile)` regrids the MxN (depth x time) data A onto an M_regularized x N grid, where the vertical postings of the new grid correspond to the `z_center` field of the table `Profile`. Interpolation is performed vertically for each timestep. 

`A_regularized = gemb_interp(Z_center,A,Profile,interp_method)` specifies a 1D interpolation method. Default is `"linear"`. 

# Example
Run a simple GEMB simulation using synthetic data and plot the resulting grid-cell center heights as blue dots. Then plot the input `Profile.z_center` values as red lines. 

```matlab
% Initialize model parameters:
ModelParam = model_initialize_parameters(output_frequency="daily");

% Generate sample data: 
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize grid:
Profile = model_initialize_column(ModelParam, ClimateForcing);

% Run GEMB: 
OutData = gemb(Profile, ClimateForcing, ModelParam);

% Get grid cell centers: 
z_center = dz2z(OutData.dz); 

% Make a 2D time grid so each z_center has a corresponding time: 
time_2D = repmat(OutData.time,size(OutData.temperature,1),1);

figure
h_gemb = plot(time_2D, z_center, '.', color='blue');
hold on
h_profile = yline(Profile.z_center, color='red');
ylabel 'Height'
title 'Grid cell center height'
legend([h_gemb(1);h_profile(1)],'gemb output','Profile.z\_center')

exportgraphics(gcf,'gemb_interp_documentation_01.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/gemb_interp_documentation_01.jpg)

Above, we see that the grid cell center locations evolve with each time step, as cells grow, compact, merge, and split in response to accumulation and densification. 

Use `gemb_interp` to interpolate the time-evolving grid onto the consistent `z_center` values of the initial `Profile`: 

```matlab
% Regrid:
temperature_gridded = gemb_interp(z_center, OutData.temperature, Profile); 

figure
pcolor(OutData.time, Profile.z_center, temperature_gridded)
shading interp
clim([250 265])
cmocean thermal % optional colormap 
cb = colorbar; 
ylabel(cb,'Temperature (K)')
box off
ylabel 'Height (m)'

exportgraphics(gcf,'gemb_interp_documentation_02.jpg',Resolution=300)
```

![](https://chadagreene.com/GEMB_figures/gemb_interp_documentation_02.jpg)

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
