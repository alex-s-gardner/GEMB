classdef test_thermo < matlab.unittest.TestCase
    
    properties
        % Grid and State
        n = 10;
        t_vec
        dz
        d
        w_surf
        re
        
        % Fluxes
        swf
        dlw
        
        % Met Data
        t_air
        v_air
        e_air
        p_air
        
        % Parameters
        dt = 3600; % 1 hour
        rho_ice = 917;
        vz = 2;
        tz = 2;
        emis = 0.98;
        ulw_delta = 0;
        emis_thresh = 10;
        emis_method = 0;
        k_method = 1;
        verbose = false;
    end
    
    methods (TestClassSetup)
        function create_mocks(tcase)
            % thermo.m depends on two other functions. To unit test thermo in isolation,
            % we create temporary mock files for these dependencies.
            
            % 1. Mock thermal_conductivity
            % Returns a constant K = 0.5 W m-1 K-1
            fid = fopen('thermal_conductivity.m', 'w');
            fprintf(fid, 'function K = thermal_conductivity(T, d, ~, ~)\n');
            fprintf(fid, '    K = 0.5 * ones(size(T));\n');
            fprintf(fid, 'end\n');
            fclose(fid);
            
            % 2. Mock turbulent_heat_flux
            % Returns zero fluxes to isolate radiative/diffusive testing
            fid = fopen('turbulent_heat_flux.m', 'w');
            fprintf(fid, 'function [shf, lhf, L] = turbulent_heat_flux(~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~)\n');
            fprintf(fid, '    shf = 0;\n');
            fprintf(fid, '    lhf = 0;\n');
            fprintf(fid, '    L = 2.5e6;\n');
            fprintf(fid, 'end\n');
            fclose(fid);
            
            % Add current folder to path to ensure mocks are found
            addpath(pwd);
            
            % Ensure we clean up later
            tcase.addTeardown(@() delete('thermal_conductivity.m'));
            tcase.addTeardown(@() delete('turbulent_heat_flux.m'));
            
            % Add source path if it exists
            try
                tcase.applyFixture(matlab.unittest.fixtures.PathFixture('../src'));
            catch
                % Fallback if running from src root
                addpath('src');
            end
        end
    end
    
    methods (TestMethodSetup)
        function setup_vectors(tcase)
            % Initialize common profile vectors
            tcase.t_vec = 260 * ones(tcase.n, 1); % Cold profile
            tcase.dz = 0.1 * ones(tcase.n, 1);    % 10cm layers
            tcase.d = 400 * ones(tcase.n, 1);     % Firn density
            tcase.w_surf = 0;
            tcase.re = 0.5 * ones(tcase.n, 1);
            
            % Initialize Fluxes (Zeroed)
            tcase.swf = zeros(tcase.n, 1);
            tcase.dlw = 0;
            
            % Met Data
            tcase.t_air = 260;
            tcase.v_air = 5;
            tcase.e_air = 100;
            tcase.p_air = 100000;
        end
    end
    
    methods (Test)
        
        function test_steady_state(tcase)
            % If T_surf = T_air and no radiation, T should remain roughly constant
            % (Assuming our mock turbulent flux is 0, and we set LW to balance)
            
            % Calculate blackbody radiation to balance LW out
            sb = 5.67e-8;
            lw_balance = sb * tcase.t_vec(1)^4 * tcase.emis;
            
            [t_out, ~, ~, ~, ~] = thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, lw_balance, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, tcase.dt, ...
                tcase.rho_ice, tcase.vz, tcase.tz, tcase.emis, tcase.ulw_delta, ...
                tcase.emis_thresh, tcase.emis_method, tcase.k_method, tcase.verbose);
            
            % Allow small diffusion drift
            tcase.verifyEqual(t_out(1), tcase.t_vec(1), 'AbsTol', 0.1, ...
                'Temperature should remain stable under balanced conditions');
        end
        
        function test_solar_heating(tcase)
            % Apply strong shortwave radiation to top cell
            tcase.swf(1) = 200; % 200 W/m2 absorbed in top layer
            
            % CRITICAL FIX: We must balance the outgoing Longwave first, 
            % otherwise the surface radiates ~250 W/m2 into space and COOLS down 
            % even with 200 W/m2 solar input.
            sb = 5.67e-8;
            lw_balance = sb * tcase.t_vec(1)^4 * tcase.emis;
            
            [t_out, shf, ~, ~, ~] = thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, lw_balance, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, tcase.dt, ...
                tcase.rho_ice, tcase.vz, tcase.tz, tcase.emis, tcase.ulw_delta, ...
                tcase.emis_thresh, tcase.emis_method, tcase.k_method, tcase.verbose);
            
            tcase.verifyTrue(t_out(1) > tcase.t_vec(1), ...
                'Top layer should warm up with SW input (given balanced LW)');
            
            % Verify internal mock works (shf should be 0)
            tcase.verifyEqual(shf, 0, 'Mocked turbulent flux should return 0');
        end
        
        function test_thermal_diffusion(tcase)
            % Set top layer hot, bottom layers cold. Turn off external fluxes.
            % Heat should diffuse to layer 2.
            
            tcase.t_vec(1) = 273;
            tcase.t_vec(2:end) = 250;
            
            % Use a long timestep to allow diffusion
            dt_long = 3600 * 3; 
            
            % CRITICAL FIX: Set emissivity to 0. 
            % If emis=0.98 and dlw=0, the hot surface radiates ~315 W/m2 loss.
            % This rapid cooling reverses the gradient before heat can diffuse down.
            emis_zero = 0;
            
            [t_out, ~, ~, ~, ~] = thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, 0, 250, 0, 0, 100000, dt_long, ...
                tcase.rho_ice, tcase.vz, tcase.tz, emis_zero, tcase.ulw_delta, ...
                tcase.emis_thresh, tcase.emis_method, tcase.k_method, tcase.verbose);
            
            % Top should cool down (giving heat to layer 2)
            tcase.verifyTrue(t_out(1) < 273, 'Hot top layer should cool via diffusion');
            
            % Layer 2 should warm up (diffusion from top)
            tcase.verifyTrue(t_out(2) > 250, 'Layer 2 should warm due to diffusion from hot surface');
        end
        
        function test_energy_conservation_check(tcase)
            % Enable verbose mode to trigger the internal energy conservation check
            % If the logic inside thermo is broken, this will throw an error.
            
            verbose_on = true;
            
            % Add some SW flux
            tcase.swf(1) = 50;
            
            try
                thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                    tcase.swf, tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, tcase.dt, ...
                    tcase.rho_ice, tcase.vz, tcase.tz, tcase.emis, tcase.ulw_delta, ...
                    tcase.emis_thresh, tcase.emis_method, tcase.k_method, verbose_on);
            catch ME
                tcase.verifyFail(['Energy conservation check failed with error: ' ME.message]);
            end
        end
        
        function test_boundary_conditions(tcase)
            % Test bottom boundary condition (Fixed T)
            % The code fixes the bottom temperature logic manually:
            % T_bottom = T(end); ... check ... if T_bottom ~= T(end) error
            
            tcase.t_vec(end) = 240; % Distinct bottom temp
            
            [t_out, ~, ~, ~, ~] = thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, tcase.dt, ...
                tcase.rho_ice, tcase.vz, tcase.tz, tcase.emis, tcase.ulw_delta, ...
                tcase.emis_thresh, tcase.emis_method, tcase.k_method, tcase.verbose);
            
            tcase.verifyEqual(t_out(end), 240, 'Bottom temperature should be fixed (Dirichlet BC)');
        end
        
        function test_timestep_subcycling(tcase)
            % thermo.m calculates a stability limit and sub-cycles if dt0 is too large.
            % We provide a very large dt0 and ensure it finishes without instability (NaNs).
            
            dt_huge = 86400; % 1 day step
            
            [t_out, ~, ~, ~, ~] = thermo(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, dt_huge, ...
                tcase.rho_ice, tcase.vz, tcase.tz, tcase.emis, tcase.ulw_delta, ...
                tcase.emis_thresh, tcase.emis_method, tcase.k_method, tcase.verbose);
            
            tcase.verifyFalse(any(isnan(t_out)), 'Solution should not explode (NaN) given large timestep');
        end
    end
end