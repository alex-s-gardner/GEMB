function [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain] = ...
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ...
    ClimateForcingStep, ModelParam, verbose)
% calculate_accumulation adds precipitation and deposition to the model column.
%
%% Syntax
%
%
%
%% Description
%
%
%
%% Inputs
%
%
%
%% Outputs
%
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% MAIN FUNCTION

% Define tolerances to allow for some numerical noise when testing equalities:  
T_tolerance    = 1e-10;
d_tolerance    = 1e-11;
gdn_tolerance  = 1e-10;
P_tolerance    = 1e-6;

% Specify constants:
CtoK              = 273.15;    % Kelvin to Celsius conversion
LF                = 0.3345E6;  % latent heat of fusion (J kg-1)
C_ice             = 2102;      % specific heat capacity of snow/ice (J kg-1 K-1)
re_new_snow       = 0.05;      % new snow grain size [mm]
gdn_new_snow      = 1.0;       % new snow dendricity
gsp_new_snow      = 0.5;       % new snow sphericity
rain              = 0;         % rainfall [mm w.e. or kg m^-2]

if verbose
    M = (dz .* density);
    M_total_initial = sum(M);                 % total mass [kg]
    E_total_initial = sum(M .* temperature * C_ice) + ...
        sum(water .* (LF + CtoK * C_ice));       % total energy [J] = initial enegy of snow/ice + initial enegy of water
end

% Density of fresh snow [kg m-3]
switch ModelParam.new_snow_method
    case "150kgm2" % Default value defined above
        density_new_snow  = 150;    % density of snow [kg m-3]

    case "350kgm2" % Density of Antarctica snow
        density_new_snow = 350.0;
        %density_new_snow = 360.0; %FirnMICE Lundin et al., 2017

    case "Fausto" % Density of Greenland snow, Fausto et al., 2018
        density_new_snow = 315.0;

        %From Vionnet et al., 2012 (Crocus)
        gdn_new_snow = min(max(1.29 - 0.17*ClimateForcingStep.wind_speed,0.20),1.0);
        gsp_new_snow = min(max(0.08*ClimateForcingStep.wind_speed + 0.38,0.5),0.9);
        re_new_snow  = max(1e-1*(gdn_new_snow/.99+(1.0-1.0*gdn_new_snow/.99).*(gsp_new_snow/.99*3.0+(1.0-gsp_new_snow/.99)*4.0))/2.0,gdn_tolerance );

    case "Kaspers" %Surface snow accumulation density from Kaspers et al., 2004, Antarctica
        %density_new_snow = alpha1 + beta1*temperature + delta1*ClimateForcingStep.precipitation_mean + epsilon1*water
        %     7.36x10-2  1.06x10-3  6.69x10-2  4.77x10-3
        density_new_snow = (7.36e-2 + 1.06e-3*min(ClimateForcingStep.temperature_air_mean,CtoK-T_tolerance ) + 6.69e-2*ClimateForcingStep.precipitation_mean/1000. + 4.77e-3*ClimateForcingStep.wind_speed_mean)*1000.;

    case "KuipersMunneke" % Kuipers Munneke and others (2015), Greenland
        density_new_snow = 481.0 + 4.834*(ClimateForcingStep.temperature_air_mean-CtoK);
end

M_surface = dz(1) * density(1);

if ClimateForcingStep.precipitation > (0 + P_tolerance)
    % Determine initial mass

    % if snow
    if ClimateForcingStep.temperature_air <= (ModelParam.rain_temperature_threshold + T_tolerance)
        z_snow = ClimateForcingStep.precipitation / density_new_snow;         % depth of snow
        dfall  = gdn_new_snow;
        sfall  = gsp_new_snow;
        refall = re_new_snow;

        % if snow depth is greater than specified min dz, new cell created
        if z_snow > ModelParam.column_dzmin+d_tolerance 

            temperature      = [ClimateForcingStep.temperature_air;    temperature]; % new cell temperature
            dz               = [z_snow                 ; dz              ]; % new cell dz
            density          = [density_new_snow       ; density         ]; % new cell density
            water            = [0                      ; water           ]; % new cell water
            albedo           = [ModelParam.albedo_snow ; albedo          ]; % new cell albedo
            albedo_diffuse   = [ModelParam.albedo_snow ; albedo_diffuse  ]; % new cell albedo_diffuse
            grain_radius     = [refall                 ; grain_radius    ]; % new cell grain size
            grain_dendricity = [dfall                  ; grain_dendricity]; % new cell grain dendricity
            grain_sphericity = [sfall                  ; grain_sphericity]; % new cell grain sphericity

            % if snow depth is less than specified minimum dz snow
        else
            M_surface_new  = M_surface + ClimateForcingStep.precipitation;                 % grid cell adjust mass
            dz(1)          = dz(1) + ClimateForcingStep.precipitation/density_new_snow;    % adjust grid cell depth
            density(1)     = M_surface_new / dz(1);                            % adjust grid cell density

            % adjust variables as a linearly weighted function of mass
            % adjust temperature (assume ClimateForcingStep.precipitation is same temp as air)
            temperature(1) = ((ClimateForcingStep.temperature_air * ClimateForcingStep.precipitation) + (temperature(1) * M_surface)) / M_surface_new;

            % adjust albedo, grain_radius, grain_dendricity & grain_sphericity
            if ModelParam.albedo_method ~= "150kgm2"
                albedo(1) = (ModelParam.albedo_snow * ClimateForcingStep.precipitation + albedo(1) * M_surface) / M_surface_new;
            end

            grain_dendricity(1) = dfall;
            grain_sphericity(1) = sfall;
            grain_radius(1)  = max(0.1*(grain_dendricity(1)/.99+(1.0-1.0*grain_dendricity(1)/.99).*(grain_sphericity(1)/.99*3.0+(1.0-grain_sphericity(1)/.99)*4.0))/2,gdn_tolerance );
        end

        % if rain
    else
        % rain is added by increasing the mass and temperature of the ice
        % of the top grid cell.  Temperatures are set artifically high to
        % account for the latent heat of fusion.  This is the same as
        % directly adding liquid water to the the snow pack surface but
        % makes the numerics easier.

        % grid cell adjust mass
        M_surface_new = M_surface + ClimateForcingStep.precipitation;

        % adjust temperature
        % liquid: must account for latent heat of fusion
        temperature(1) = ((ClimateForcingStep.precipitation * (ClimateForcingStep.temperature_air + LF/C_ice)) + ...
            (temperature(1) * M_surface)) / M_surface_new;

        % adjust grid cell density
        density(1) = M_surface_new / dz(1);

        % if density > the density of ice, density = ModelParam.density_ice
        if density(1) > ModelParam.density_ice-d_tolerance 
           density(1)  = ModelParam.density_ice;          % adjust density
           dz(1) = M_surface_new / density(1);            % dz is adjusted to conserve mass
        end

        rain = ClimateForcingStep.precipitation;
    end

    if verbose
        % Check for conservation of mass:
        M             = (dz .* density);
        M_total_final = sum(M);                   % total mass [kg]
        M_delta       = M_total_final - M_total_initial - ClimateForcingStep.precipitation;

        E_total_final = sum(M .* temperature * C_ice) + ...
            sum(water .* (LF + CtoK * C_ice));       % total energy [J] = initial enegy of snow/ice + initial enegy of water

        E_snow = ((ClimateForcingStep.precipitation - rain) * (ClimateForcingStep.temperature_air ) * C_ice);
        E_rain = (rain * (ClimateForcingStep.temperature_air * C_ice + LF));

        E_delta = E_total_final - E_total_initial - E_snow - E_rain;
        
        if (abs(M_delta) > 1E-3) || (abs(E_delta) > 1E-3)
            error(['Mass and/or energy are not conserved:' newline ' M_delta: ' ...
                num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
        end
    end

end