# `gemb_spinup` documentation
`gemb_spinup` spins up a `gemb` `Profile`.

See also [`model_initialize_profile`](model_initialize_profile_documentation.md) and [`gemb_profile `](gemb_profile_documentation.md).

# Syntax 
```matlab
 Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, spinup_cycles)
 Profile_spunup = gemb_spinup(..., display_waitbar=false)
```

# Description
`Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, spinup_cycles)` runs the gemb function an integer number of spinup_cycles. Default `spinup_cycles` = 1. 

`Profile_spunup = gemb_spinup(..., display_waitbar=false)` disables the graphical waitbar. 

# Example

```matlab
% Generate sample data: 
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize model parameters:
ModelParam = model_initialize_parameters(output_frequency="daily");
   
% Initialize grid:
Profile_initial = model_initialize_profile(ModelParam, ClimateForcing);
  
% Convert full forcing time series to climatology:
ClimateForcingSpinup = forcing_climatology(ClimateForcing); 

% Spinup a profile for 50 climatological average years: 
Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, 50);

% Run GEMB: 
OutData = gemb(Profile_spunup, ClimateForcing, ModelParam);
```

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 
