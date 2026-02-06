classdef test_turbulent_heat_flux < matlab.unittest.TestCase
    
    properties
        % Standard atmospheric conditions
        temperature_air_std     = 270;    % Kelvin
        temperature_surface_std = 265;    % Kelvin (Stable condition temperature_air > T_surf typically)
        pressure_air_std        = 100000; % Pa
        vapor_pressure_std      = 400;    % Pa (Vapor pressure)
        wind_speed_std          = 5;      % m/s
        rho_air_std             = 1.25;   % kg/m3
        
        % Heights and roughness
        wind_observation_height        = 2;  % Wind measurement height
        temperature_observation_height = 2;  % Temp measurement height
        z0 = 0.001;           % Roughness
        zt = 0.0001;
        zq = 0.0001;
        
        % Physical constants
        lv = 2.495E6;
        ls = 2.8295e6;
        
        % Structs
        CF % ClimateForcingStep
    end
    
    methods (TestMethodSetup)
        function setup_structs(tcase)
            % Initialize ClimateForcingStep (CF) with standard values
            tcase.CF.temperature_air                = tcase.temperature_air_std;
            tcase.CF.pressure_air                   = tcase.pressure_air_std;
            tcase.CF.vapor_pressure                 = tcase.vapor_pressure_std;
            tcase.CF.wind_speed                     = tcase.wind_speed_std;
            tcase.CF.wind_observation_height        = tcase.wind_observation_height;
            tcase.CF.temperature_observation_height = tcase.temperature_observation_height;
            
            % Add source path
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
        end
    end
    
    methods (Test)
        
        function test_stable_conditions(tcase)
            % Stable: temperature_air > T_surface (Ri > 0)
            tcase.CF.temperature_air = 275;
            temperature_surface      = 270; % Cold surface
            
            [heat_flux_sensible, heat_flux_latent, ~] = turbulent_heat_flux(temperature_surface, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            % In stable conditions with temperature_air > T_surf, sensible heat flows INTO surface (+SHF)
            tcase.verifyTrue(heat_flux_sensible > 0, 'SHF should be positive (downward) for warm air over cold surface');
            
            % Verify coefficients didn't explode
            tcase.verifyFalse(isnan(heat_flux_sensible));
            tcase.verifyFalse(isnan(heat_flux_latent));
        end
        
        function test_unstable_conditions(tcase)
            % Unstable: T_surface > temperature_air (Ri < 0)
            tcase.CF.temperature_air = 260;
            temperature_surface      = 270; % Hot surface
            
            [heat_flux_sensible, heat_flux_latent, ~] = turbulent_heat_flux(temperature_surface, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            % In unstable conditions with T_surf > temperature_air, sensible heat flows AWAY from surface (-SHF)
            tcase.verifyTrue(heat_flux_sensible < 0, 'SHF should be negative (upward) for cold air over warm surface');
        end
        
        function test_sublimation_ice_phase(tcase)
            % T_surface < 273.15 -> Ice -> Sublimation (L = LS)
            temperature_surface = 260; 
            
            [~, ~, l_out] = turbulent_heat_flux(temperature_surface, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            % Verify correct latent heat constant
            tcase.verifyEqual(l_out, tcase.ls, 'AbsTol', 1e-1, 'Should use Latent Heat of Sublimation for ice');
        end
        
        function test_evaporation_water_phase(tcase)
            % T_surface >= 273.15 -> Water -> Vaporization (L = LV)
            temperature_surface = 275; 
            
            [~, ~, l_out] = turbulent_heat_flux(temperature_surface, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            % Verify correct latent heat constant
            tcase.verifyEqual(l_out, tcase.lv, 'AbsTol', 1e-1, 'Should use Latent Heat of Vaporization for water');
        end
        
        function test_latent_heat_direction(tcase)
            % Test dry air over wet surface -> Evaporation -> Negative LHF (Heat loss)
            tcase.CF.vapor_pressure = 0; % Very dry air
            temperature_surface     = 273.15; % Saturation pressure ~611 Pa
            
            [~, heat_flux_latent, ~] = turbulent_heat_flux(temperature_surface, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            tcase.verifyTrue(heat_flux_latent < 0, 'Evaporation should result in negative latent heat flux (surface cooling)');
            
            % Test humid air over dry/cold surface -> Condensation -> Positive LHF (Heat gain)
            tcase.CF.vapor_pressure = 1000; 
            temperature_surface_cold = 250; % Saturation pressure very low
            
            [~, heat_flux_latent, ~] = turbulent_heat_flux(temperature_surface_cold, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            tcase.verifyTrue(heat_flux_latent > 0, 'Condensation/Deposition should result in positive latent heat flux');
        end
        
        function test_zero_wind_stability(tcase)
            % Wind speed 0 usually causes singularities in Monin-Obukhov. 
            % Verify the function handles small/zero wind without crashing (NaN/Inf).
            
            tcase.CF.wind_speed = 0.0001; % Close to zero
            
            [heat_flux_sensible, heat_flux_latent, ~] = turbulent_heat_flux(tcase.temperature_surface_std, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            tcase.verifyFalse(isnan(heat_flux_sensible) || isinf(heat_flux_sensible), 'Fluxes should not be NaN/Inf at near-zero wind');
            tcase.verifyTrue(abs(heat_flux_sensible) < 1.0, 'Fluxes should be negligible at zero wind');
        end
        
        function test_neutral_stability_check(tcase)
            % Neutral: temperature_air == T_surface -> Ri ~ 0
            % In neutral conditions, coefM should approach log(z/z0).
            
            t_iso = 270;
            tcase.CF.temperature_air = t_iso;
            
            [heat_flux_sensible, ~, ~] = turbulent_heat_flux(t_iso, tcase.rho_air_std, ...
                tcase.z0, tcase.zt, tcase.zq, tcase.CF);
            
            tcase.verifyEqual(heat_flux_sensible, 0, 'AbsTol', 1e-10, 'Sensible heat flux must be 0 if dT is 0');
        end
    end
end