%% estimate simulation parameters from real data 
% simulation_parameters_estimate_from_data Calculates statistical coefficients for synthetic climate generation.
%
%% Syntax
%
% simulation_parameters_estimate_from_data
%
%% Description
%
% This script analyzes real-world meteorological data to derive the statistical 
% coefficients required by the GEMB synthetic climate generator. It performs the 
% following operations:
%
% 1. Data Ingestion: Loads a reference dataset (e.g., TEST_INPUT_1.mat) 
%    containing observed time series for temperature, radiation, wind, etc..
% 2. Parameter Fitting: Fits statistical models to the observed data to 
%    extract key generative parameters [Image of statistical parameter estimation]:
%    * Temperature: Fits mean, seasonal amplitude, and daily noise characteristics.
%    * Humidity & Wind: Fits autoregressive noise models to residuals after 
%      removing seasonal trends.
%    * Precipitation: Fits harmonic functions to occurrence probabilities 
%      and magnitude distributions.
%    * Longwave Radiation: Fits a Gaussian mixture model to the residuals 
%      between observed longwave and clear-sky estimates (representing cloud effects).
% 3. Validation: Simulates synthetic data using the derived coefficients and 
%    plots it alongside the original observations for visual verification.
% 4. Output Generation: Prints the derived coefficients to the console in 
%    executable MATLAB format, ready to be pasted into simulation_parameter_sets.m.
%
%% Inputs
%
%  fn (hardcoded) : string       Path to the input .mat file containing observed climate data structure inputs.
%
%% Outputs
%
%  Console Output : text         Formatted MATLAB code block defining location_parameters and coeffs.
%  Figures        : plots        Comparison plots of Observed vs. Simulated time series for each variable.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% estimate simulation parameters from real data 
fn = '/Users/gardnera/Code/GEMB/GEMB_0.21/TEST/TEST_INPUT_1.mat'; % path to input data 
inputs = load(fn);

%% location and time parameters
location_parameters.description = "parameters estimated using simulation_parameters_estimate_from_data.m as fit to original TEST_INPUT_1.mat data";

location_parameters.latitude = inputs.LP.latitude;
location_parameters.longitude = inputs.LP.longitude;
location_parameters.elevation = 700; %meters 
location_parameters.time_step = 1/(365.25*24); % fraction of a year
location_parameters.start_date = 1994; % decimal year
location_parameters.end_date = location_parameters.start_date + 31; % decimal year
location_parameters.rand_seed = 42; % Sets the seed to a fixed number

location_parameters.wind_observation_height = inputs.LP.wind_observation_height; % wind observation height above surface [m]
location_parameters.temperature_observation_height = inputs.LP.temperature_observation_height; % temperature observation height above surface [m]
location_parameters.temperature_air_mean = mean(inputs.temperature_air); % average annual temerature [K]
location_parameters.precipitation_mean = mean(inputs.precipitation); % average annual accumulation rate of snow or ice [kg m⁻² yr⁻¹]

dec_year = location_parameters.start_date: location_parameters.time_step:location_parameters.end_date+1;
dec_year = dec_year(1:length(inputs.temperature_air));
dec_year = dec_year(:);
rng(location_parameters.rand_seed)

disp("%%  location and time parameters")
disp("location_parameters.description = """ + location_parameters.description + """;")
disp("location_parameters.latitude = " + sprintf('%0.4f ', location_parameters.latitude)+ "; % [º]")
disp("location_parameters.longitude = " + sprintf('%0.4f ', location_parameters.longitude)+ "; % [º]")
disp("location_parameters.elevation = " + sprintf('%0.0f ', location_parameters.elevation)+ "; % [m]")
disp("location_parameters.start_date = " + sprintf('%0.2f ', location_parameters.start_date) + "; % [decimal year]")
disp("location_parameters.end_date = " + sprintf('%0.2f ', location_parameters.end_date) + "; % [decimal year]")

disp("location_parameters.wind_observation_height = " + sprintf('%0.1f ', location_parameters.wind_observation_height)+ "; % wind observation height above surface [m]")
disp("location_parameters.temperature_observation_height = " + sprintf('%0.1f ', location_parameters.temperature_observation_height)+ "; % temperature observation height above surface [m]")
disp("location_parameters.temperature_air_mean = " + sprintf('%0.1f ', location_parameters.temperature_air_mean)+ "; % average annual temerature [K]")
disp("location_parameters.precipitation_mean = " + sprintf('%0.1f ', location_parameters.precipitation_mean)+ "; % average annual accumulation rate of snow or ice [kg m⁻² yr⁻¹]")

disp("location_parameters.time_step = " + sprintf('%0.4f ', location_parameters.time_step) + "; % [fraction of a year]")
disp("location_parameters.rand_seed = " + sprintf('%0.0f ', 42) + "; % [seed for random number generator]")
disp(" ")

%% downward shortwave [W/m2]
varname = "shortwave_downward";
longname = varname2longname(varname);
simulated.(varname) = simulate_shortwave_irradiance(dec_year, location_parameters.latitude);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level air temperature [K]
varname = "temperature_air";
longname = varname2longname(varname);
disp("%% " + longname)
coeffs.(varname) = fit_air_temperature(dec_year, inputs.(varname + "0"), location_parameters.latitude, location_parameters.elevation);
simulated.(varname ) = simulate_air_temperature(dec_year, location_parameters.latitude, location_parameters.elevation, coeffs.(varname));
simulate_coeffs_disp(coeffs.(varname), "coeffs." +varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level air pressure [Pa]
varname = "pressure_air";
longname = varname2longname(varname);
simulated.(varname)= simulate_air_pressure(dec_year, simulated.temperature_air, location_parameters.latitude, location_parameters.elevation);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level relative humidity
varname = "relative_humidity";
min_max = [0, 100]';
longname = varname2longname(varname);
disp("%% " + longname)
inputs.(varname + "0") = relative_humidity(inputs.vapor_pressure, inputs.temperature_air);
coeffs.(varname) = fit_seasonal_daily_noise(dec_year, inputs.(varname + "0"));
coeffs.(varname).min_max = min_max;
simulate_coeffs_disp(coeffs.(varname), "coeffs." + varname)
simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
simulated.(varname)(simulated.(varname)<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
simulated.(varname)(simulated.(varname)>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
simulate_coeffs_disp(coeffs.(varname), "coeffs." +varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level vapor pressure [Pa]
varname = "vapor_pressure";
longname = varname2longname(varname);
simulated.(varname) = simulate_vapor_pressure(simulated.temperature_air, simulated.relative_humidity);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% downward longwave [W/m2]
varname = "longwave_downward";
min_max = [0, Inf]';
longname = varname2longname(varname);
disp("%% " + longname)
simulated.(varname) = simulate_longwave_irradiance(simulated.temperature_air, simulated.vapor_pressure);

% account for cloud cover 
coeffs.(varname) = fit_longwave_irradiance_delta(inputs.(varname + "0") - simulated.(varname));
simulated.(varname) = simulated.(varname)  + simulate_longwave_irradiance_delta(dec_year, coeffs.(varname));
coeffs.(varname).min_max = min_max;
simulated.(varname)(simulated.(varname)<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
simulated.(varname)(simulated.(varname)>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
simulate_coeffs_disp(coeffs.(varname), "coeffs." + varname)

figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname));plot(inputs.(varname + "0")); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen wind speed [m/s]
varname = "wind_speed";
min_max = [0, Inf]';
longname = varname2longname(varname);
disp("%% " + longname)
coeffs.(varname) = fit_seasonal_daily_noise(dec_year, inputs.(varname + "0"));
coeffs.(varname).min_max = min_max;
simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
simulated.(varname)(simulated.(varname)<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
simulated.(varname)(simulated.(varname)>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
simulate_coeffs_disp(coeffs.(varname), "coeffs."+varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% Precipitation [kg m-2]
varname = "precipitation";
min_max = [0, Inf]';
longname = varname2longname(varname);
disp("%% " + longname)
coeffs.(varname) = fit_precipitation(dec_year, inputs.(varname + "0"));
simulated.(varname) = simulate_precipitation(dec_year, coeffs.(varname));
simulate_coeffs_disp(coeffs.(varname), "coeffs."+varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); plot(inputs.(varname + "0"));
ylabel(longname); legend(["observed", "simulated"]); hold off;