function [ClimateForcing] = simulate_climate_forcing(set_id)

    % load climate simulation parameter set
    [location_parameters, coeffs] = simulation_parameter_sets(set_id);
    
    % initialize times and random seed
    dec_year = location_parameters.start_date:location_parameters.time_step:location_parameters.end_date+1;
    dec_year = dec_year(:);
    daten = decyear2datenum(dec_year);
    rng(location_parameters.rand_seed);
    
    % simulate downward shortave radiation
    dsw0 = simulate_shortwave_irradiance(dec_year, location_parameters.lat);
    
    % simulate downward longwave radiation
    T_air0 = simulate_air_temperature(dec_year, location_parameters.lat, location_parameters.elev, coeffs.T_air);
    
    % screen level air temperature [K]
    p_air0 = simulate_air_pressure(dec_year, T_air0, location_parameters.lat, location_parameters.elev);
    
    % screen level relative humidity [%]
    varname = "rh";
    rh0 = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
    rh0(rh0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
    rh0(rh0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
    
    % screen level vapor pressure [Pa]
    e_air0 = simulate_vapor_pressure(T_air0, rh0);
    
    % downward logwave radiation [W m⁻²]
    varname = "dlw";
    dlw0 = simulate_longwave_irradiance(T_air0, e_air0);
    dlw0 = dlw0  + simulate_longwave_irradiance_delta(dec_year, coeffs.(varname));
    dlw0(dlw0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
    dlw0(dlw0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
    
    % screen level wind speed [m s⁻¹]
    varname = "V";
    V0 = simulate_seasonal_daily_noise(dec_year, coeffs.(varname));
    V0(V0<coeffs.(varname).min_max(1)) = coeffs.(varname).min_max(1);
    V0(V0>coeffs.(varname).min_max(2)) = coeffs.(varname).min_max(2);
    

    % precipitation [kg m⁻²]
    varname = "P";
    P0 = simulate_precipitation(dec_year, coeffs.(varname));

    % populate structure
    ClimateForcing.daten = daten;
    ClimateForcing.dsw0 = dsw0;
    ClimateForcing.T_air0 = T_air0;
    ClimateForcing.p_air0 = p_air0;
    ClimateForcing.rh0 = rh0;
    ClimateForcing.e_air0 = e_air0;
    ClimateForcing.dlw0 = dlw0;
    ClimateForcing.V0 = V0;
    ClimateForcing.P0 = P0;

    % Location specifc parameters
    ClimateForcing.Vz = location_parameters.Vz;
    ClimateForcing.Tz = location_parameters.Tz;
    ClimateForcing.T_air_mean = location_parameters.T_air_mean;
    ClimateForcing.V_mean = mean(ClimateForcing.V0);
    ClimateForcing.P_mean = location_parameters.P_mean;
    ClimateForcing.elev = location_parameters.elev;
    ClimateForcing.lat = location_parameters.lat;
    ClimateForcing.lon = location_parameters.lon;
end