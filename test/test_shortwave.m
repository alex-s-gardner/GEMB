classdef test_shortwave < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 5;
        dz
        d
        re
        dsw = 200;       % Downward SW [W m-2]
        dsw_diff = 50;   % Diffuse SW [W m-2]
        alb = 0.7;       % Surface albedo
        alb_diff = 0.8;  % Diffuse albedo
        rho_ice = 917;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize profiles
            tcase.dz = 0.1 * ones(tcase.n, 1);
            tcase.d = 350 * ones(tcase.n, 1); % Snow density
            tcase.re = 0.5 * ones(tcase.n, 1); % 0.5 mm grain radius
            
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
            % Albedo method != 1 (Standard calculation)
            
            method_sw = 0;
            method_alb = 3; 
            
            swf = shortwave(tcase.dz, tcase.d, tcase.re, tcase.dsw, tcase.dsw_diff, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, method_alb, method_sw);
            
            % Expect all energy in top cell
            expected_net = (1 - tcase.alb) * tcase.dsw;
            
            tcase.verifyEqual(swf(1), expected_net, 'AbsTol', 1e-6, ...
                'Top cell should absorb all net SW when absorption method is 0');
            tcase.verifyEqual(sum(swf(2:end)), 0, 'AbsTol', 1e-6, ...
                'Lower cells should receive 0 energy');
        end
        
        function test_surface_absorption_gardner(tcase)
            % Test Case: No penetration, Gardner Albedo Method (method_alb = 1)
            % This uses a specific formula separating diffuse and direct components
            
            method_sw = 0;
            method_alb = 1;
            
            swf = shortwave(tcase.dz, tcase.d, tcase.re, tcase.dsw, tcase.dsw_diff, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, method_alb, method_sw);
            
            dsw_direct = tcase.dsw - tcase.dsw_diff;
            expected = (1 - tcase.alb) * dsw_direct + (1 - tcase.alb_diff) * tcase.dsw_diff;
            
            tcase.verifyEqual(swf(1), expected, 'AbsTol', 1e-6, ...
                'Gardner method (1) should treat diffuse and direct components separately');
        end
        
        function test_density_override(tcase)
            % Test Case: Penetration requested (sw_method = 1), but surface is ICE
            % Should revert to surface absorption
            
            method_sw = 1;
            method_alb = 3;
            
            d_ice = tcase.d;
            d_ice(1) = tcase.rho_ice; % Set top cell to ice density
            
            swf = shortwave(tcase.dz, d_ice, tcase.re, tcase.dsw, tcase.dsw_diff, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, method_alb, method_sw);
            
            expected = (1 - tcase.alb) * tcase.dsw;
            
            tcase.verifyEqual(swf(1), expected, 'AbsTol', 1e-6, ...
                'Should revert to surface absorption if density >= ice density');
            tcase.verifyEqual(sum(swf(2:end)), 0, 'AbsTol', 1e-6);
        end
        
        function test_penetration_brun_method(tcase)
            % Test Case: Penetration with Brun method (method_alb = 2)
            % This method calculates its own spectral albedos internally
            
            method_sw = 1;
            method_alb = 2;
            
            swf = shortwave(tcase.dz, tcase.d, tcase.re, tcase.dsw, tcase.dsw_diff, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, method_alb, method_sw);
            
            % 1. Energy must be distributed deeper than top cell
            tcase.verifyTrue(swf(2) > 0, 'Brun method should penetrate below surface');
            
            % 2. Profile shape: Absorption should generally decrease with depth 
            % (assuming uniform dz and density)
            tcase.verifyTrue(swf(1) > swf(2), 'Absorption should decrease with depth');
            tcase.verifyTrue(swf(2) > swf(3));
            
            % 3. Total absorption should be less than incoming (due to albedo)
            tcase.verifyTrue(sum(swf) < tcase.dsw, 'Total absorbed cannot exceed downward flux');
            
            % 4. Conservation check (internal consistency)
            tcase.verifyTrue(sum(swf) > 0.05 * tcase.dsw, 'Should absorb a significant fraction');
        end
        
        function test_penetration_standard_method(tcase)
            % Test Case: Penetration with Standard method (method_alb != 2)
            % Uses Greuell & Konzelmann coefficients
            
            method_sw = 1;
            method_alb = 3;
            
            % CRITICAL FIX: Use a deep column to prevent flux from escaping the bottom.
            % Original n=5 (0.5m) allowed ~0.6% energy to exit.
            % Increasing to 100 layers (10m) ensures ~100% absorption.
            n_deep = 100;
            dz_deep = 0.1 * ones(n_deep, 1);
            d_deep = 350 * ones(n_deep, 1);
            re_deep = 0.5 * ones(n_deep, 1);
            
            swf = shortwave(dz_deep, d_deep, re_deep, tcase.dsw, tcase.dsw_diff, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, method_alb, method_sw);
            
            % 1. Conservation: Sum of absorbed flux must equal Surface Net Flux
            expected_total = (1 - tcase.alb) * tcase.dsw;
            
            tcase.verifyEqual(sum(swf), expected_total, 'AbsTol', 1e-4, ...
                'Total absorbed energy must be conserved (using deep column to capture all flux)');
            
            % 2. Distribution Check
            % Standard method: 36% absorbed at surface (band > 0.8um) + 
            % penetrating part absorbed at surface.
            tcase.verifyTrue(swf(1) > 0.36 * expected_total, ...
                'Top cell should absorb NIR band + penetrating UV/Vis part');
            
            % 3. Penetration Check
            tcase.verifyTrue(swf(2) > 0, 'Energy should penetrate to second layer');
        end
        
        function test_zero_flux(tcase)
            % Test Case: Night time (0 incoming SW)
            
            swf = shortwave(tcase.dz, tcase.d, tcase.re, 0, 0, ...
                tcase.alb, tcase.alb_diff, tcase.rho_ice, 3, 1);
            
            tcase.verifyEqual(sum(swf), 0, 'AbsTol', 1e-10, 'Zero input flux should result in zero absorption');
        end
    end
end