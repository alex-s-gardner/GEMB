function [T_air, V, dlw, dsw, e_air, p_air, P, black_carbon_snow, ...
    black_carbon_ice, cloud_optical_thickness, solar_zenith_angle, ...
    dsw_diffuse, cloud_fraction] = ...
   model_inputs_single_timestep(index, T_air0, V0, dlw0, dsw0, e_air0, ...
   p_air0, P0, S)
        

        T_air   = T_air0(index);    % screen level air temperature [K]   
        V       = V0(index);        % wind speed [m s-1]
        dlw     = dlw0(index);      % downward longwave radiation flux [W m-2]
        dsw     = dsw0(index);      % downward shortwave radiation flux [W m-2]
        e_air   = e_air0(index);    % screen level vapor pressure [Pa]
        p_air   = p_air0(index);    % screen level air pressure [Pa]
        P       = P0(index);        % precipitation [kg m-2]
      
       

        % if we are provided with cc and cot values, extract for the timestep
        if numel(S.black_carbon_snow)>1
            black_carbon_snow = S.black_carbon_snow(index);
        else
            black_carbon_snow = S.black_carbon_snow;
        end

        if numel(S.black_carbon_ice)>1
            black_carbon_ice = S.black_carbon_ice(index);
        else
            black_carbon_ice = S.black_carbon_ice;
        end

        if numel(S.cloud_optical_thickness)>1
            cloud_optical_thickness = S.cloud_optical_thickness(index);
        else
            cloud_optical_thickness = S.cloud_optical_thickness;
        end

        if numel(S.solar_zenith_angle)>1
            solar_zenith_angle = S.solar_zenith_angle(index);
        else
            solar_zenith_angle = S.solar_zenith_angle;
        end

        if numel(S.dsw_diffuse)>1
            dsw_diffuse = S.dsw_diffuse(index);
        else
            dsw_diffuse = S.dsw_diffuse;
        end

        if numel(S.cloud_fraction)>1
            cloud_fraction = S.cloud_fraction(index);
        else
            cloud_fraction = S.cloud_fraction;
        end
