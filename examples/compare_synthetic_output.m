% Compare MATLAB synthetic example output with Julia version
% This script runs the synthetic example and outputs comparable statistics

addpath("../src")

% Generate 3-hourly synthetic climate forcing data:
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize model parameters:
ModelParam = model_initialize_parameters(output_frequency="daily");

% Initialize a column:
Profile = model_initialize_profile(ModelParam, ClimateForcing);

% Create a climatological average time series:
ClimateForcingClimatology = forcing_climatology(ClimateForcing);

% Spinup a profile for 75 years of average forcing:
ModelParam_spinup = model_initialize_parameters(output_frequency="last");
Profile_spunup = gemb_spinup(Profile, ClimateForcingClimatology, ModelParam_spinup, 75);

% Run GEMB with the spun-up profile:
OutData = gemb(Profile_spunup, ClimateForcing, ModelParam);

% Print summary statistics matching Julia output:
fprintf('Simulation complete!\n');
fprintf('  Time steps: %d\n', length(OutData.time));
fprintf('  Profile layers: %d\n', size(OutData.temperature, 1));
fprintf('  Mean surface albedo: %.3f\n', mean(OutData.albedo_surface, 'omitnan'));
fprintf('  Total melt: %.2f kg/m²\n', sum(OutData.melt, 'omitnan'));
fprintf('  Total runoff: %.2f kg/m²\n', sum(OutData.runoff, 'omitnan'));

% Save key outputs for detailed comparison
save('matlab_synthetic_output.mat', 'OutData', '-v7.3');
fprintf('\nOutput saved to matlab_synthetic_output.mat\n');
