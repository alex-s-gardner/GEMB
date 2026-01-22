classdef test_shortwave < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 5;
        dz
        d
        re
        
        alb = 0.7;       % Surface albedo
        alb_diff = 0.8;  % Diffuse albedo
        
        % Structs
        CF % ClimateForcingStep
        MP % ModelParam
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize profiles
            tcase.dz = 0.1 * ones(tcase.n, 1);
            tcase.d = 350 * ones(tcase.n, 1); % Snow density
            tcase.re = 0.5 * ones(tcase.n, 1); % 0.5 mm grain radius
            
            % Initialize ClimateForcingStep (CF)
            tcase.CF.dsw = 200;       % Downward SW [W m-2]
            tcase.CF.dsw_diffuse = 50;   % Diffuse SW [W m-2]
            
            % Initialize ModelParam (MP)
            tcase.MP.density_ice = 917;
            tcase.MP.sw_absorption_method = 0; % Default: No penetration
            tcase.MP.albedo_method = "GreuellKonzelmann"; % Default
            
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
        
        function test_surface_absorption_basic(tcase)
            % Test Case: No penetration (sw_absorption_method = 0)
            % Albedo method != "GardnerSharp" (Standard calculation)
            
            tcase.MP.sw_absorption_method = 0;
            tcase.MP.albedo_method = "GreuellKonzelmann"; 
            
            swf = calculate_shortwave_radiation(tcase.dz, tcase.d, tcase.re, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            % Expect all energy in top cell
            expected_net = (1 - tcase.alb) * tcase.CF.dsw;
            
            tcase.verifyEqual(swf(1), expected_net, 'AbsTol', 1e-6, ...
                'Top cell should absorb all net SW when absorption method is 0');
            tcase.verifyEqual(sum(swf(2:end)), 0, 'AbsTol', 1e-6, ...
                'Lower cells should receive 0 energy');
        end
        
        function test_surface_absorption_gardner(tcase)
            % Test Case: No penetration, Gardner Albedo Method ("GardnerSharp")
            % This uses a specific formula separating diffuse and direct components
            
            tcase.MP.sw_absorption_method = 0;
            tcase.MP.albedo_method = "GardnerSharp";
            
            swf = calculate_shortwave_radiation(tcase.dz, tcase.d, tcase.re, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            dsw_direct = tcase.CF.dsw - tcase.CF.dsw_diffuse;
            expected = (1 - tcase.alb) * dsw_direct + (1 - tcase.alb_diff) * tcase.CF.dsw_diffuse;
            
            tcase.verifyEqual(swf(1), expected, 'AbsTol', 1e-6, ...
                'Gardner method should treat diffuse and direct components separately');
        end
        
        function test_density_override(tcase)
            % Test Case: Penetration requested (sw_method = 1), but surface is ICE
            % Should revert to surface absorption
            
            tcase.MP.sw_absorption_method = 1;
            tcase.MP.albedo_method = "GreuellKonzelmann";
            
            d_ice = tcase.d;
            d_ice(1) = tcase.MP.density_ice; % Set top cell to ice density
            
            swf = calculate_shortwave_radiation(tcase.dz, d_ice, tcase.re, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            expected = (1 - tcase.alb) * tcase.CF.dsw;
            
            tcase.verifyEqual(swf(1), expected, 'AbsTol', 1e-6, ...
                'Should revert to surface absorption if density >= ice density');
            tcase.verifyEqual(sum(swf(2:end)), 0, 'AbsTol', 1e-6);
        end
        
        function test_penetration_brun_method(tcase)
            % Test Case: Penetration with Brun method ("BruneLeFebre")
            % This method calculates its own spectral albedos internally
            
            tcase.MP.sw_absorption_method = 1;
            tcase.MP.albedo_method = "BruneLeFebre";
            
            swf = calculate_shortwave_radiation(tcase.dz, tcase.d, tcase.re, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            % 1. Energy must be distributed deeper than top cell
            tcase.verifyTrue(swf(2) > 0, 'Brun method should penetrate below surface');
            
            % 2. Profile shape: Absorption should generally decrease with depth 
            tcase.verifyTrue(swf(1) > swf(2), 'Absorption should decrease with depth');
            tcase.verifyTrue(swf(2) > swf(3));
            
            % 3. Total absorption should be less than incoming (due to albedo)
            tcase.verifyTrue(sum(swf) < tcase.CF.dsw, 'Total absorbed cannot exceed downward flux');
            
            % 4. Conservation check (internal consistency)
            tcase.verifyTrue(sum(swf) > 0.05 * tcase.CF.dsw, 'Should absorb a significant fraction');
        end
        
        function test_penetration_standard_method(tcase)
            % Test Case: Penetration with Standard method (NOT "BruneLeFebre")
            % Uses Greuell & Konzelmann coefficients
            
            tcase.MP.sw_absorption_method = 1;
            tcase.MP.albedo_method = "GreuellKonzelmann";
            
            % CRITICAL FIX: Use a deep column to prevent flux from escaping the bottom.
            n_deep = 100;
            dz_deep = 0.1 * ones(n_deep, 1);
            d_deep = 350 * ones(n_deep, 1);
            re_deep = 0.5 * ones(n_deep, 1);
            
            swf = calculate_shortwave_radiation(dz_deep, d_deep, re_deep, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            % 1. Conservation: Sum of absorbed flux must equal Surface Net Flux
            expected_total = (1 - tcase.alb) * tcase.CF.dsw;
            
            tcase.verifyEqual(sum(swf), expected_total, 'AbsTol', 1e-4, ...
                'Total absorbed energy must be conserved (using deep column to capture all flux)');
            
            % 2. Distribution Check
            tcase.verifyTrue(swf(1) > 0.36 * expected_total, ...
                'Top cell should absorb NIR band + penetrating UV/Vis part');
            
            % 3. Penetration Check
            tcase.verifyTrue(swf(2) > 0, 'Energy should penetrate to second layer');
        end
        
        function test_zero_flux(tcase)
            % Test Case: Night time (0 incoming SW)
            tcase.CF.dsw = 0;
            tcase.CF.dsw_diffuse = 0;
            tcase.MP.sw_absorption_method = 1;
            
            swf = calculate_shortwave_radiation(tcase.dz, tcase.d, tcase.re, tcase.alb, tcase.alb_diff, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyEqual(sum(swf), 0, 'AbsTol', 1e-10, 'Zero input flux should result in zero absorption');
        end
    end
end