function [T, dz, d, water, re, gdn, gsp, a, a_diffuse, rain] = ...
    calculate_accumulation(T, dz, d, water, re, gdn, gsp, a, a_diffuse, ...
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
CI                = 2102;      % specific heat capacity of snow/ice (J kg-1 k-1)
re_new_snow       = 0.05;      % new snow grain size [mm]
gdn_new_snow      = 1.0;       % new snow dendricity
gsp_new_snow      = 0.5;       % new snow sphericity
rain              = 0;         % rainfall [mm w.e. or kg m^-3]

if verbose
    M = (dz .* d);
    M_total_initial = sum(M);                 % total mass [kg]
    E_total_initial = sum(M .* T * CI) + ...
        sum(water .* (LF + CtoK * CI));       % total energy [J] = initial enegy of snow/ice + initial enegy of water
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
        gdn_new_snow = min(max(1.29 - 0.17*ClimateForcingStep.V,0.20),1.0);
        gsp_new_snow = min(max(0.08*ClimateForcingStep.V + 0.38,0.5),0.9);
        re_new_snow  = max(1e-1*(gdn_new_snow/.99+(1.0-1.0*gdn_new_snow/.99).*(gsp_new_snow/.99*3.0+(1.0-gsp_new_snow/.99)*4.0))/2.0,gdn_tolerance );

    case "Kaspers" %Surface snow accumulation density from Kaspers et al., 2004, Antarctica
        %density_new_snow = alpha1 + beta1*T + delta1*ClimateForcingStep.P_mean + epsilon1*water
        %     7.36x10-2  1.06x10-3  6.69x10-2  4.77x10-3
        density_new_snow = (7.36e-2 + 1.06e-3*min(ClimateForcingStep.T_air_mean,CtoK-T_tolerance ) + 6.69e-2*ClimateForcingStep.P_mean/1000. + 4.77e-3*ClimateForcingStep.V_mean)*1000.;

    case "KuipersMunneke" % Kuipers Munneke and others (2015), Greenland
        density_new_snow = 481.0 + 4.834*(ClimateForcingStep.T_air_mean-CtoK);
end

M_surface = dz(1) * d(1);

if ClimateForcingStep.P > (0 + P_tolerance)
    % Determine initial mass

    % if snow
    if ClimateForcingStep.T_air <= (ModelParam.rain_temperature_threshold + CtoK + T_tolerance)
        z_snow = ClimateForcingStep.P / density_new_snow;         % depth of snow
        dfall  = gdn_new_snow;
        sfall  = gsp_new_snow;
        refall = re_new_snow;

        % if snow depth is greater than specified min dz, new cell created
        if z_snow > ModelParam.column_dzmin+d_tolerance 

            T         = [ClimateForcingStep.T_air;    T        ]; % new cell T
            dz        = [z_snow                  ;    dz       ]; % new cell dz
            d         = [density_new_snow        ;    d        ]; % new cell d
            water     = [0                       ;    water    ]; % new cell water
            a         = [ModelParam.albedo_snow  ;    a        ]; % new cell a
            a_diffuse = [ModelParam.albedo_snow  ;    a_diffuse]; % new cell a_diffuse
            re        = [refall                  ;    re       ]; % new cell grain size
            gdn       = [dfall                   ;    gdn      ]; % new cell grain dendricity
            gsp       = [sfall                   ;    gsp      ]; % new cell grain sphericity

            % if snow depth is less than specified minimum dz snow
        else
            M_surface_new  = M_surface + ClimateForcingStep.P;                 % grid cell adjust mass
            dz(1)          = dz(1) + ClimateForcingStep.P/density_new_snow;    % adjust grid cell depth
            d(1)           = M_surface_new / dz(1);                            % adjust grid cell density

            % adjust variables as a linearly weighted function of mass
            % adjust temperature (assume ClimateForcingStep.P is same temp as air)
            T(1) = ((ClimateForcingStep.T_air * ClimateForcingStep.P) + (T(1) * M_surface)) / M_surface_new;

            % adjust a, re, gdn & gsp
            if ModelParam.albedo_method ~= "150kgm2"
                a(1) = (ModelParam.albedo_snow * ClimateForcingStep.P + a(1) * M_surface) / M_surface_new;
            end

            gdn(1) = dfall;
            gsp(1) = sfall;
            re(1)  = max(0.1*(gdn(1)/.99+(1.0-1.0*gdn(1)/.99).*(gsp(1)/.99*3.0+(1.0-gsp(1)/.99)*4.0))/2,gdn_tolerance );
        end

        % if rain
    else
        % rain is added by increasing the mass and temperature of the ice
        % of the top grid cell.  Temperatures are set artifically high to
        % account for the latent heat of fusion.  This is the same as
        % directly adding liquid water to the the snow pack surface but
        % makes the numerics easier.

        % grid cell adjust mass
        M_surface_new = M_surface + ClimateForcingStep.P;

        % adjust temperature
        % liquid: must account for latent heat of fusion
        T(1) = ((ClimateForcingStep.P * (ClimateForcingStep.T_air + LF/CI)) + ...
            (T(1) * M_surface)) / M_surface_new;

        % adjust grid cell density
        d(1) = M_surface_new / dz(1);

        % if d > the density of ice, d = ModelParam.density_ice
        if d(1) > ModelParam.density_ice-d_tolerance 
           d(1)  = ModelParam.density_ice;          % adjust d
           dz(1) = M_surface_new / d(1);            % dz is adjusted to conserve mass
        end

        rain = ClimateForcingStep.P;
    end

    if verbose
        % Check for conservation of mass:
        M = (dz .* d);
        M_total_final = sum(M);                   % total mass [kg]
        M_delta = M_total_final - M_total_initial - ClimateForcingStep.P;

        E_total_final = sum(M .* T * CI) + ...
            sum(water .* (LF + CtoK * CI));       % total energy [J] = initial enegy of snow/ice + initial enegy of water

        E_snow = ((ClimateForcingStep.P - rain) * (ClimateForcingStep.T_air ) * CI);
        E_rain = (rain * (ClimateForcingStep.T_air * CI + LF));

        E_delta = E_total_final - E_total_initial - E_snow - E_rain;
        
        if (abs(M_delta) > 1E-3) || (abs(E_delta) > 1E-3)
            error(['Mass and/or energy are not conserved:' newline ' M_delta: ' ...
                num2str(M_delta) ' E_delta: ' num2str(E_delta) newline])
        end
    end

end