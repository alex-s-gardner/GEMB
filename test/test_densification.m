classdef test_densification < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 10;
        t_vec
        dz
        d
        re
        
        % Structs
        CF % ClimateForcingStep
        MP % ModelParam
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
            
            % Initialize ClimateForcingStep (CF)
            tcase.CF.dt = 86400 * 30; % 30 days in seconds
            tcase.CF.P_mean = 200;    % kg m-2 yr-1
            tcase.CF.T_air_mean = 250; % Kelvin
            
            % Initialize ModelParam (MP)
            tcase.MP.density_ice = 917;
            tcase.MP.densification_method = "HerronLangway"; % Default
            tcase.MP.densification_coeffs_M01 = "Ant_RACMO_GS_SW0"; % Default
            
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
            
            tcase.MP.densification_method = "HerronLangway";
            
            [dz_out, d_out] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
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
            d_near_ice(:) = tcase.MP.density_ice - 0.1; 
            
            % Long timestep to force overshoot
            tcase.CF.dt = 86400 * 365 * 100; 
            
            tcase.MP.densification_method = "HerronLangway";
            
            [~, d_out] = calculate_density(tcase.t_vec, tcase.dz, d_near_ice, tcase.re, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(d_out <= tcase.MP.density_ice + 1e-10), ...
                'Density must be clamped at density_ice');
        end
        
        function test_herron_langway(tcase)
            % Test "HerronLangway" logic
            tcase.MP.densification_method = "HerronLangway";
            [~, d_out] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'HerronLangway should densify layers');
        end
        
        function test_arthern(tcase)
            % Test "Anthern" (Semi-empirical Arthern)
            tcase.MP.densification_method = "Anthern";
            [~, d_out] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'Anthern should densify layers');
        end
        
        function test_arthern_b_physical(tcase)
            % "AnthernB" uses grain radius (re) and overburden pressure
            tcase.MP.densification_method = "AnthernB";
            
            % Run standard
            [~, d_std] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            % Run with larger grain size
            re_large = tcase.re * 2;
            [~, d_large] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, re_large, ...
                tcase.CF, tcase.MP);
            
            % Formula B1 in Arthern: rate is proportional to 1/r^2. 
            % Larger grains -> Slower densification -> Lower final density
            % Note: Check applies where d < density_ice
            
            diff = d_std - d_large;
            
            % Exclude top layer (overburden = 0) and clamped values
            valid_mask = (d_std < tcase.MP.density_ice - 1);
            
            % Create an index list to check, strictly excluding index 1
            check_indices = find(valid_mask);
            check_indices(check_indices == 1) = [];
            
            if ~isempty(check_indices)
                tcase.verifyTrue(all(diff(check_indices) > 0), ...
                    'AnthernB: Larger grains should densify slower (excluding surface layer)');
            end
        end
        
        function test_ligtenberg_coeffs(tcase)
            % Test "Ligtenberg" method with different coefficients
            tcase.MP.densification_method = "Ligtenberg";
            
            % Test Case 1: Standard RACMO
            tcase.MP.densification_coeffs_M01 = "Ant_RACMO_GS_SW0";
            [~, d1] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            % Test Case 2: ERA5 variant (should have different M coeffs)
            tcase.MP.densification_coeffs_M01 = "Ant_ERA5_BF_SW1";
            [~, d2] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            % Ensure densification occurred
            tcase.verifyTrue(all(d1 > tcase.d));
            tcase.verifyTrue(all(d2 > tcase.d));
            
            % Verify different coeffs produce different results
            % (Assuming the M lookup values differ for these keys)
            tcase.verifyFalse(isequal(d1, d2), 'Different Ligtenberg M01 coeffs should yield different results');
        end
        
        function test_kuipers_munneke_coeffs(tcase)
            % Test Greenland logic via Ligtenberg method with specific coeffs
            tcase.MP.densification_method = "Ligtenberg";
            tcase.MP.densification_coeffs_M01 = "Gre_KuipersMunneke";
            
            [~, d_out] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(d_out > tcase.d), 'Gre_KuipersMunneke coeffs should densify layers');
        end
        
        function test_zero_time_step(tcase)
            % If dt = 0, no change should occur
            tcase.CF.dt = 0;
            tcase.MP.densification_method = "HerronLangway";
            
            [dz_out, d_out] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            tcase.verifyEqual(d_out, tcase.d);
            tcase.verifyEqual(dz_out, tcase.dz);
        end
        
        function test_ligtenberg_bare_ice_logic(tcase)
            % Test the specific logic branch in "Ligtenberg" for 820 vs 917 density_ice
            tcase.MP.densification_method = "Ligtenberg";
            tcase.MP.densification_coeffs_M01 = "Gre_RACMO_GS_SW0"; % Requires 2-row M01 matrix
            
            % Case A: density_ice ~ 820 (Trigger specialized branch)
            tcase.MP.density_ice = 820;
            [~, d_820] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            % Case B: density_ice ~ 917 (Trigger standard branch)
            tcase.MP.density_ice = 917;
            [~, d_917] = calculate_density(tcase.t_vec, tcase.dz, tcase.d, tcase.re, ...
                tcase.CF, tcase.MP);
            
            % The code uses M01(1,:) for 820 and M01(2,:) for others.
            % Assuming Gre_RACMO_GS_SW0 returns a 2x4 matrix where rows differ.
            
            % Normalize result for comparison (since d_max clamping differs)
            % Check only values well below 820
            check_idx = tcase.d < 700;
            
            tcase.verifyFalse(isequal(d_820(check_idx), d_917(check_idx)), ...
                'Density_ice=820 logic branch should produce different rates than standard');
        end
    end
end