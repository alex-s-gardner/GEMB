
%% estimate simulation parameters from real data 
fn = '/Users/gardnera/Code/GEMB/GEMB_0.21/TEST/TEST_INPUT_1.mat'; % path to input data 
inputs = load(fn);

%% location and time parameters
location_parameters.lat = -73.3307;
location_parameters.elev = 700; %meters 
location_parameters.time_step = 1/(365.25*24); % fraction of a year
location_parameters.start_date = 1994; % decimal year
location_parameters.end_date = location_parameters.start_date + 31; % decimal year
location_parameters.rand_seed = 42; % Sets the seed to a fixed number

dec_year = location_parameters.start_date: location_parameters.time_step:location_parameters.end_date+1;
dec_year = dec_year(1:length(inputs.Ta0));
rand(location_parameters.rand_seed);

disp("%%  location and time parameters")
disp("location_parameters.lat = " + sprintf('%0.4f ', location_parameters.lat)+ "; % [º]")
disp("location_parameters.elev = " + sprintf('%0.0f ', location_parameters.elev)+ "; % [m]")
disp("location_parameters.start_date = " + sprintf('%0.2f ', location_parameters.start_date) + "; % [decimal year]")
disp("location_parameters.end_date = " + sprintf('%0.2f ', location_parameters.end_date) + "; % [decimal year]")
disp("location_parameters.time_step = " + sprintf('%0.4f ', location_parameters.time_step) + "; % [fraction of a year]")
disp("location_parameters.rand_seed = " + sprintf('%0.0f ', 42) + "; % [seed for random number generator]")
disp(" ")

%% downward shortwave [W/m2]
varname = "dsw";
longname = varname2longname(varname);
simulated.(varname) = simulate_shortwave_irradiance(dec_year, location_parameters.lat);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level air temperature [K]
varname = "Ta";
longname = varname2longname(varname);
disp("%% " + longname)
coeffs.(varname) = fit_air_temperature(dec_year, inputs.(varname + "0"), location_parameters.lat, location_parameters.elev);
simulated.(varname ) = simulate_air_temperature(dec_year, location_parameters.lat, location_parameters.elev, coeffs.(varname));
simulate_coeffs_disp(coeffs.(varname), "coeffs.(varname)." +varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level air pressure [Pa]
varname = "pAir";
longname = varname2longname(varname);
simulated.(varname)= simulate_air_pressure(dec_year,simulated.Ta, location_parameters.lat, location_parameters.elev);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level relative humidity
varname = "rh";
min_max = [0, 100]';
longname = varname2longname(varname);
disp("%% " + longname)
inputs.(varname + "0") = relative_humidity(inputs.eAir0, inputs.Ta0);
coeffs.(varname) = fit_seasonal_daily_noise(dec_year, inputs.(varname + "0"));
coeffs.(varname).min_max = min_max;
simulate_coeffs_disp(coeffs.(varname), "coeffs.(varname)." + varname)
simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
simulated.(varname) (simulated.(varname)<min_max(1)) = min_max(1);
simulated.(varname) (simulated.(varname)>min_max(2)) = min_max(2);
simulate_coeffs_disp(coeffs.(varname), "coeffs.(varname)." +varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen level vapor pressure [Pa]
varname = "eAir";
longname = varname2longname(varname);
simulated.(varname) = simulate_vapor_pressure(simulated.Ta, simulated.rh);
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% downward longwave [W/m2]
varname = "dlw";
longname = varname2longname(varname);
disp("%% " + longname)
simulated.(varname) = simulate_longwave_irradiance(simulated.Ta, simulated.eAir);

% account for cloud cover 
coeffs.(varname) = fit_seasonal_daily_noise(dec_year, inputs.(varname + "0") - simulated.(varname));
simulated.(varname) = simulated.(varname)  + simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
simulate_coeffs_disp(coeffs.(varname), "coeffs.(varname)." + varname)

figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;

%% screen wind speed [m/s]
varname = "V";
min_max = [0, Inf]';
longname = varname2longname(varname);
disp("%% " + longname)
coeffs.(varname) = fit_seasonal_daily_noise(dec_year, inputs.(varname + "0"));
coeffs.(varname).min_max = min_max;
simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
simulated.(varname) (simulated.(varname)<min_max(1)) = min_max(1);
simulated.(varname) (simulated.(varname)>min_max(2)) = min_max(2);
simulate_coeffs_disp(coeffs.(varname), "coeffs.(varname)." +varname)
figure; plot(inputs.(varname + "0")); hold on; plot(simulated.(varname)); 
ylabel(longname); legend(["observed", "simulated"]); hold off;