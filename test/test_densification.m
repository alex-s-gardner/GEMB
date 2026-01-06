classdef test_densification < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 10;
        t_vec
        dz
        d
        re
        dt = 86400 * 30; % 30 days in seconds
        
        % Constants
        rho_ice = 917;
        t_mean = 250;
        p_mean = 200; % kg m-2 yr-1
        
        % Calibration flags (Defaults)
        alb_dens_thresh = 1023;
        alb_method = 1;
        sw_method = 0;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize profile
            % T: Cold to prevent melt, but allow densification
            tcase.t_vec = 260 * ones(tcase.n, 1);
            tcase.dz = 0.5 * ones(tcase.n, 1);
            
            % D: Mixed profile crossing the 550 kg/m3 threshold
            % Top half < 550, Bottom half > 550
            tcase.d = linspace(300, 700, tcase.n)';
            
            tcase.re = 0.5 * ones(tcase.n, 1);
            
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
        
        function test_mass_conservation(tcase)
            % Densification changes density (d) and thickness (dz).
            % Mass (d * dz) must remain constant.
            
            mass_initial = tcase.d .* tcase.dz;
            
            method = 1; % Herron and Langway
            
            [dz_out, d_out] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            mass_final = d_out .* dz_out;
            
            tcase.verifyEqual(mass_final, mass_initial, 'AbsTol', 1e-10, ...
                'Mass must be conserved during densification');
            
            % Check that densification actually happened
            tcase.verifyTrue(all(d_out > tcase.d), 'Density should increase over time');
            tcase.verifyTrue(all(dz_out < tcase.dz), 'Thickness should decrease over time');
        end
        
        function test_density_clamping(tcase)
            % Density should never exceed ice density
            
            % Set density very close to ice density
            d_near_ice = tcase.d;
            d_near_ice(:) = tcase.rho_ice - 0.1; 
            
            % Long timestep to force overshoot
            dt_long = 86400 * 365 * 100; 
            
            method = 1;
            
            [~, d_out] = densification(tcase.t_vec, tcase.dz, d_near_ice, tcase.re, ...
                dt_long, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            tcase.verifyTrue(all(d_out <= tcase.rho_ice + 1e-10), ...
                'Density must be clamped at density_ice');
        end
        
        function test_herron_langway_method_1(tcase)
            % Test Method 1 logic
            method = 1;
            [~, d_out] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'Method 1 should densify layers');
        end
        
        function test_arthern_method_2(tcase)
            % Test Method 2 (Semi-empirical Arthern)
            method = 2;
            [~, d_out] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'Method 2 should densify layers');
        end
        
        function test_arthern_physical_method_3(tcase)
            % Method 3 uses grain radius (re) and overburden pressure
            method = 3;
            
            % Run standard
            [~, d_std] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            % Run with larger grain size
            re_large = tcase.re * 2;
            [~, d_large] = densification(tcase.t_vec, tcase.dz, tcase.d, re_large, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            % Formula B1 in Arthern: rate is proportional to 1/r^2. 
            % Larger grains -> Slower densification -> Lower final density
            % Note: Check applies where d < density_ice
            
            diff = d_std - d_large;
            
            % Filter out clamped values
            % CRITICAL FIX: Also exclude the top layer (index 1).
            % Method 3 depends on overburden pressure, which is 0 for the top layer.
            % Therefore, densification is 0 for the top layer regardless of grain size.
            valid_mask = (d_std < tcase.rho_ice - 1);
            
            % Create an index list to check, strictly excluding index 1
            check_indices = find(valid_mask);
            check_indices(check_indices == 1) = [];
            
            if ~isempty(check_indices)
                tcase.verifyTrue(all(diff(check_indices) > 0), ...
                    'Method 3: Larger grains should densify slower (excluding surface layer)');
            end
        end
        
        function test_ligtenberg_calibration_method_6(tcase)
            % Method 6 has specific branches for ERA5 calibration
            method = 6;
            
            % Branch 1: albedo=1, sw=0 (Standard)
            [~, d1] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, 1, 0, method);
            
            % Branch 2: albedo=2, sw=1
            [~, d2] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, 2, 1, method);
            
            % Branch 3: Default (RACMO)
            [~, d3] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, 2, 0, method);
            
            % Ensure all branches produced valid densification different from input
            tcase.verifyTrue(all(d1 > tcase.d));
            tcase.verifyTrue(all(d2 > tcase.d));
            tcase.verifyTrue(all(d3 > tcase.d));
            
            % Verify branches produce distinct results (coefficients differ)
            tcase.verifyFalse(isequal(d1, d2), 'Different calibration branches should yield different results');
        end
        
        function test_kuipers_munneke_method_7(tcase)
            % Method 7 (Greenland)
            method = 7;
            
            [~, d_out] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'Method 7 should densify layers');
        end
        
        function test_zero_time_step(tcase)
            % If dt = 0, no change should occur
            dt_zero = 0;
            method = 1;
            
            [dz_out, d_out] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                dt_zero, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                tcase.alb_dens_thresh, tcase.alb_method, tcase.sw_method, method);
            
            tcase.verifyEqual(d_out, tcase.d);
            tcase.verifyEqual(dz_out, tcase.dz);
        end
        
        function test_bare_ice_calibration_branch(tcase)
            % Test the specific sub-branch in Method 6 for "bare ice" calibration
            % Triggered when abs(albedo_desnity_threshold - 820.0) >= tolerance
            method = 6;
            
            % Case A: Near 820 (Standard)
            alb_thresh_std = 820;
            [~, d_std] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                alb_thresh_std, 1, 0, method);
            
            % Case B: Far from 820 (Bare Ice Logic)
            alb_thresh_ice = 1023;
            [~, d_ice] = densification(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.dt, tcase.rho_ice, tcase.t_mean, tcase.p_mean, ...
                alb_thresh_ice, 1, 0, method);
            
            % Coefficients differ, so results should differ
            tcase.verifyFalse(isequal(d_std, d_ice), 'Bare ice calibration branch should differ from standard');
        end
    end
end