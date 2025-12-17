function simulated = simulate_gemb_forcing(location_parameters, coeffs)
    %% simulate_gemb_forcing 
    
    dec_year = location_parameters.start_date: location_parameters.time_step:location_parameters.end_date+1;
    dec_year = dec_year(1:length(inputs.Ta0));
    rand(location_parameters.rand_seed);
    
    %% downward shortwave [W/m2]
    varname = "dsw";
    simulated.(varname) = simulate_shortwave_irradiance(dec_year, location_parameters.lat);
    
    %% screen level air temperature [K]
    varname = "Ta";
    simulated.(varname ) = simulate_air_temperature(dec_year, location_parameters.lat, location_parameters.elev, coeffs.(varname));
    
    %% screen level air pressure [Pa]
    varname = "pAir";
    simulated.(varname)= simulate_air_pressure(dec_year,simulated.Ta, location_parameters.lat, location_parameters.elev);
    
    %% screen level relative humidity
    varname = "rh";
    simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
    simulated.(varname) (simulated.(varname)<min_max(1)) = min_max(1);
    simulated.(varname) (simulated.(varname)>min_max(2)) = min_max(2);
    
    %% screen level vapor pressure [Pa]
    varname = "eAir";
    simulated.(varname) = simulate_vapor_pressure(simulated.Ta, simulated.rh);
    
    %% downward longwave [W/m2]
    varname = "dlw";
    simulated.(varname) = simulate_longwave_irradiance(simulated.Ta, simulated.eAir);
    
    % account for cloud cover 
    simulated.(varname) = simulated.(varname)  + simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
    
    %% screen wind speed [m/s]
    varname = "V";
    simulated.(varname) = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
    simulated.(varname) (simulated.(varname)<min_max(1)) = min_max(1);
    simulated.(varname) (simulated.(varname)>min_max(2)) = min_max(2);
end
