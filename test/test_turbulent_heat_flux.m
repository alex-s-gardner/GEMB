classdef test_turbulent_heat_flux < matlab.unittest.TestCase
    
    properties
        % Standard atmospheric conditions
        t_air_std = 270;      % Kelvin
        t_surf_std = 265;     % Kelvin (Stable condition T_air > T_surf typically)
        p_air_std = 100000;   % Pa
        e_air_std = 400;      % Pa (Vapor pressure)
        v_std = 5;            % m/s
        rho_air_std = 1.25;   % kg/m3
        
        % Heights and roughness
        vz = 2;               % Wind measurement height
        tz = 2;               % Temp measurement height
        z0 = 0.001;           % Roughness
        zt = 0.0001;
        zq = 0.0001;
        
        % Physical constants
        lv = 2.495E6;
        ls = 2.8295e6;
    end
    
    methods (TestMethodSetup)
        function add_source_path(tcase)
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
            % Stable: T_air > T_surface (Ri > 0)
            t_air = 275;
            t_surf = 270; % Cold surface
            
            [shf, lhf, ~] = turbulent_heat_flux(t_air, t_surf, tcase.p_air_std, ...
                tcase.e_air_std, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            % In stable conditions with T_air > T_surf, sensible heat flows INTO surface (+SHF)
            tcase.verifyTrue(shf > 0, 'SHF should be positive (downward) for warm air over cold surface');
            
            % Verify coefficients didn't explode
            tcase.verifyFalse(isnan(shf));
            tcase.verifyFalse(isnan(lhf));
        end
        
        function test_unstable_conditions(tcase)
            % Unstable: T_surface > T_air (Ri < 0)
            t_air = 260;
            t_surf = 270; % Hot surface
            
            [shf, lhf, ~] = turbulent_heat_flux(t_air, t_surf, tcase.p_air_std, ...
                tcase.e_air_std, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            % In unstable conditions with T_surf > T_air, sensible heat flows AWAY from surface (-SHF)
            tcase.verifyTrue(shf < 0, 'SHF should be negative (upward) for cold air over warm surface');
        end
        
        function test_sublimation_ice_phase(tcase)
            % T_surface < 273.15 -> Ice -> Sublimation (L = LS)
            t_surf = 260; 
            
            [~, ~, l_out] = turbulent_heat_flux(tcase.t_air_std, t_surf, tcase.p_air_std, ...
                tcase.e_air_std, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            % Verify correct latent heat constant
            tcase.verifyEqual(l_out, tcase.ls, 'AbsTol', 1e-1, 'Should use Latent Heat of Sublimation for ice');
        end
        
        function test_evaporation_water_phase(tcase)
            % T_surface >= 273.15 -> Water -> Vaporization (L = LV)
            t_surf = 275; 
            
            [~, ~, l_out] = turbulent_heat_flux(tcase.t_air_std, t_surf, tcase.p_air_std, ...
                tcase.e_air_std, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            % Verify correct latent heat constant
            tcase.verifyEqual(l_out, tcase.lv, 'AbsTol', 1e-1, 'Should use Latent Heat of Vaporization for water');
        end
        
        function test_latent_heat_direction(tcase)
            % Test dry air over wet surface -> Evaporation -> Negative LHF (Heat loss)
            e_air_dry = 0; % Very dry air
            t_surf = 273.15; % Saturation pressure ~611 Pa
            
            [~, lhf, ~] = turbulent_heat_flux(tcase.t_air_std, t_surf, tcase.p_air_std, ...
                e_air_dry, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            tcase.verifyTrue(lhf < 0, 'Evaporation should result in negative latent heat flux (surface cooling)');
            
            % Test humid air over dry/cold surface -> Condensation -> Positive LHF (Heat gain)
            e_air_humid = 1000; 
            t_surf_cold = 250; % Saturation pressure very low
            
            [~, lhf, ~] = turbulent_heat_flux(tcase.t_air_std, t_surf_cold, tcase.p_air_std, ...
                e_air_humid, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            tcase.verifyTrue(lhf > 0, 'Condensation/Deposition should result in positive latent heat flux');
        end
        
        function test_zero_wind_stability(tcase)
            % Wind speed 0 usually causes singularities in Monin-Obukhov. 
            % Verify the function handles small/zero wind without crashing (NaN/Inf).
            % Note: Calling code usually clamps V, but if passed directly:
            
            v_zero = 0.0001; % Close to zero
            
            [shf, lhf, ~] = turbulent_heat_flux(tcase.t_air_std, tcase.t_surf_std, tcase.p_air_std, ...
                tcase.e_air_std, v_zero, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            tcase.verifyFalse(isnan(shf) || isinf(shf), 'Fluxes should not be NaN/Inf at near-zero wind');
            tcase.verifyTrue(abs(shf) < 1.0, 'Fluxes should be negligible at zero wind');
        end
        
        function test_neutral_stability_check(tcase)
            % Neutral: T_air == T_surface -> Ri ~ 0
            % In neutral conditions, coefM should approach log(z/z0).
            
            t_iso = 270;
            
            % We cannot check the internal variables (coefM) directly as they are not returned 
            % in the standard output list of the provided function signature 
            % [shf, lhf, L]. 
            % However, we can verify SHF is 0 when T_air == T_surf.
            
            [shf, ~, ~] = turbulent_heat_flux(t_iso, t_iso, tcase.p_air_std, ...
                tcase.e_air_std, tcase.v_std, tcase.rho_air_std, tcase.vz, tcase.tz, ...
                tcase.z0, tcase.zt, tcase.zq);
            
            tcase.verifyEqual(shf, 0, 'AbsTol', 1e-10, 'Sensible heat flux must be 0 if dT is 0');
        end
    end
end