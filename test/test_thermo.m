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
        
        % Structures
        CF % ClimateForcingStep
        MP % ModelParam
        
        verbose = false;
    end
    
    methods (TestClassSetup)
        function create_mocks(tcase)
            % calculate_temperature.m depends on helper functions. To unit test calculate_temperature in isolation,
            % we create temporary mock files for these dependencies.
            
            % 1. Mock thermal_conductivity
            % Returns a constant K = 0.5 W m-1 K-1
            fid = fopen('thermal_conductivity.m', 'w');
            fprintf(fid, 'function K = thermal_conductivity(T, d, ~)\n');
            fprintf(fid, '    K = 0.5 * ones(size(T));\n');
            fprintf(fid, 'end\n');
            fclose(fid);
            
            % 2. Mock turbulent_heat_flux
            % Returns zero fluxes to isolate radiative/diffusive testing
            fid = fopen('turbulent_heat_flux.m', 'w');
            fprintf(fid, 'function [shf, lhf, L] = turbulent_heat_flux(~, ~, ~, ~, ~, ~)\n');
            fprintf(fid, '    shf = 0;\n');
            fprintf(fid, '    lhf = 0;\n');
            fprintf(fid, '    L = 2.5e6;\n');
            fprintf(fid, 'end\n');
            fclose(fid);

            % 3. Mock thermo_optimal_dt
            % FIX: Return a stable timestep (3600s). 
            % Previous version returned max(divisors) = 86400, which caused 
            % numerical instability (checking Von Neumann stability) and test failure.
            fid = fopen('thermo_optimal_dt.m', 'w');
            fprintf(fid, 'function dt = thermo_optimal_dt(~, ~, ~, ~, ~)\n');
            fprintf(fid, '    dt = 3600;\n'); 
            fprintf(fid, 'end\n');
            fclose(fid);
            
            % Add current folder to path to ensure mocks are found
            addpath(pwd);
            
            % Ensure we clean up later
            tcase.addTeardown(@() delete('thermal_conductivity.m'));
            tcase.addTeardown(@() delete('turbulent_heat_flux.m'));
            tcase.addTeardown(@() delete('thermo_optimal_dt.m'));
            
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
        function setup_inputs(tcase)
            % Initialize common profile vectors
            tcase.t_vec = 260 * ones(tcase.n, 1); % Cold profile
            tcase.dz = 0.1 * ones(tcase.n, 1);    % 10cm layers
            tcase.d = 400 * ones(tcase.n, 1);     % Firn density
            tcase.w_surf = 0;
            tcase.re = 0.5 * ones(tcase.n, 1);
            
            % Initialize Fluxes (Zeroed)
            tcase.swf = zeros(tcase.n, 1);
            
            % Initialize ClimateForcingStep (CF)
            tcase.CF.dt = 3600; % 1 hour
            tcase.CF.dlw = 0;
            tcase.CF.T_air = 260;
            tcase.CF.V = 5;
            tcase.CF.e_air = 100;
            tcase.CF.p_air = 100000;
            tcase.CF.Vz = 2;
            tcase.CF.Tz = 2;
            
            % Initialize ModelParam (MP)
            tcase.MP.density_ice = 917;
            tcase.MP.emissivity = 0.98;
            tcase.MP.emissivity_re_large = 0.98; % Default fallback
            tcase.MP.ulw_delta = 0;
            tcase.MP.emissivity_re_threshold = 10;

            tcase.MP.surface_roughness_effective_ratio = 0.1; 
            
            % calculate_temperature.m expects a string for the method
            tcase.MP.emissivity_method = "uniform"; 
            
            tcase.MP.thermal_conductivity_method = "Sturm"; % Mocked anyway
            
            % calculate_temperature.m requires dt_divisors to determine stability
            tcase.MP.dt_divisors = [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, ...
                                    120, 300, 600, 900, 1200, 1800, 3600, 86400];
        end
    end
    
    methods (Test)
        
        function test_steady_state(tcase)
            % If T_surf = T_air and no radiation, T should remain roughly constant
            % (Assuming our mock turbulent flux is 0, and we set LW to balance)
            
            % Calculate blackbody radiation to balance LW out
            sb = 5.67e-8;
            tcase.CF.dlw = sb * tcase.t_vec(1)^4 * tcase.MP.emissivity;
            
            [t_out, ~, ~, ~, ~] = calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.CF, tcase.MP, tcase.verbose);
            
            % Allow small diffusion drift
            tcase.verifyEqual(t_out(1), tcase.t_vec(1), 'AbsTol', 0.1, ...
                'Temperature should remain stable under balanced conditions');
        end
        
        function test_solar_heating(tcase)
            % Apply strong shortwave radiation to top cell
            tcase.swf(1) = 200; % 200 W/m2 absorbed in top layer
            
            % Balance the outgoing Longwave first so surface doesn't cool
            sb = 5.67e-8;
            tcase.CF.dlw = sb * tcase.t_vec(1)^4 * tcase.MP.emissivity;
            
            [t_out, ~, shf, ~, ~] = calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.CF, tcase.MP, tcase.verbose);
            
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
            tcase.CF.dt = 3600 * 3; 
            
            % Turn off emissivity to prevent radiative cooling at surface
            tcase.MP.emissivity = 0;
            tcase.CF.dlw = 0;
            tcase.CF.T_air = 250;
            
            [t_out, ~, ~, ~, ~] = calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.CF, tcase.MP, tcase.verbose);
            
            % Top should cool down (giving heat to layer 2)
            tcase.verifyTrue(t_out(1) < 273, 'Hot top layer should cool via diffusion');
            
            % Layer 2 should warm up (diffusion from top)
            tcase.verifyTrue(t_out(2) > 250, 'Layer 2 should warm due to diffusion from hot surface');
        end
        
        function test_energy_conservation_check(tcase)
            % Enable verbose mode to trigger the internal energy conservation check
            verbose_on = true;
            
            % Add some SW flux
            tcase.swf(1) = 50;
            
            try
                calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                    tcase.swf, tcase.CF, tcase.MP, verbose_on);
            catch ME
                tcase.verifyFail(['Energy conservation check failed with error: ' ME.message]);
            end
        end
        
        function test_boundary_conditions(tcase)
            % Test bottom boundary condition (Fixed T)
            tcase.t_vec(end) = 240; % Distinct bottom temp
            
            [t_out, ~, ~, ~, ~] = calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(t_out(end), 240, 'Bottom temperature should be fixed (Dirichlet BC)');
        end
        
        function test_timestep_subcycling(tcase)
            % calculate_temperature.m calculates a stability limit and sub-cycles if dt is too large.
            % We provide a very large dt (1 day) and ensure it finishes without instability (NaNs).
            % The mock thermo_optimal_dt returns 3600, so this will subcycle 24 times.
            
            tcase.CF.dt = 86400; % 1 day step
            
            [t_out, ~, ~, ~, ~] = calculate_temperature(tcase.t_vec, tcase.dz, tcase.d, tcase.w_surf, tcase.re, ...
                tcase.swf, tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyFalse(any(isnan(t_out)), 'Solution should not explode (NaN) given large timestep');
        end
    end
end