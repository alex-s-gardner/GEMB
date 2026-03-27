% This is a simple example of how to run GEMB using synthetic forcing. 


%% Set up the model and run it: 
addpath("src")

% Generate 3-hourly synthetic climate forcing data: 
time_step_hours = 3;
ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);

% Initialize model parameters:
ModelParam = model_initialize_parameters();

% Initialize a column:
Profile = model_initialize_profile(ModelParam, ClimateForcing);

% Spinup model (i.e. allow profile state to equilibrate to climate)
ModelParam.output_frequency = "last";
spinup_cycles = 3;
for i = 1:spinup_cycles
    OutData = gemb(Profile, ClimateForcing, ModelParam);
    Profile = gemb_profile(OutData);
end

% Run GEMB using equilibrated Profile (Takes a minute):
Profile = gemb_profile(OutData);
OutData = gemb(Profile, ClimateForcing, ModelParam);

%% Plot results: 

% Get a 2D matrix of grid cell centers: 
z_center = dz2z(OutData.dz);

% Convert time to 2D so pcolor can plot it: 
time_2D = repmat(OutData.time,size(OutData.temperature,1),1);

figure
pcolor(time_2D,z_center,OutData.temperature)
shading interp
clim([250 270])
ylabel 'Column height (m)'
ylim([-10 1])
cb = colorbar;
ylabel(cb,'Temperature (K)')