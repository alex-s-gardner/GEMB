classdef test_thermal_conductivity < matlab.unittest.TestCase
    
    properties
        rho_ice = 917;
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
        
        function test_sturm_method_snow(tcase)
            % Test Method 1 (Sturm et al., 1997) for snow densities
            % Formula: 0.138 - 1.01E-3*d + 3.233E-6*d^2
            
            d_snow = 300;
            t_in = 260; % Temperature shouldn't affect snow formula in this model
            method = 1;
            
            k_out = thermal_conductivity(t_in, d_snow, tcase.rho_ice, method);
            
            expected = 0.138 - 1.01e-3 * d_snow + 3.233e-6 * d_snow^2;
            
            tcase.verifyEqual(k_out, expected, 'AbsTol', 1e-8, ...
                'Method 1 should follow Sturm parameterization for snow');
        end
        
        function test_calonne_method_snow(tcase)
            % Test Method 2 (Calonne et al., 2011) for snow densities
            % Formula: 0.024 - 1.23E-4*d + 2.5e-6*d^2
            
            d_snow = 300;
            t_in = 260;
            method = 2;
            
            k_out = thermal_conductivity(t_in, d_snow, tcase.rho_ice, method);
            
            expected = 0.024 - 1.23e-4 * d_snow + 2.5e-6 * d_snow^2;
            
            tcase.verifyEqual(k_out, expected, 'AbsTol', 1e-8, ...
                'Method 2 should follow Calonne parameterization for snow');
        end
        
        function test_ice_conductivity(tcase)
            % Test Ice Physics (Density >= Density_Ice)
            % Formula: 9.828 * exp(-5.7E-3 * T)
            % Should depend on T, not density (as long as d >= threshold)
            
            d_ice = 917;
            t_cold = 240;
            t_warm = 270;
            method = 1; % Method flag shouldn't change ice physics
            
            % Test Cold Ice
            k_cold = thermal_conductivity(t_cold, d_ice, tcase.rho_ice, method);
            expected_cold = 9.828 * exp(-5.7e-3 * t_cold);
            tcase.verifyEqual(k_cold, expected_cold, 'AbsTol', 1e-8);
            
            % Test Warm Ice
            k_warm = thermal_conductivity(t_warm, d_ice, tcase.rho_ice, method);
            expected_warm = 9.828 * exp(-5.7e-3 * t_warm);
            tcase.verifyEqual(k_warm, expected_warm, 'AbsTol', 1e-8);
            
            % Verify T dependence
            tcase.verifyNotEqual(k_cold, k_warm, 'Ice conductivity should depend on temperature');
        end
        
        function test_mixed_profile(tcase)
            % Test a vector containing both snow and ice
            t_vec = [260; 250];
            d_vec = [400; 920]; % 400=Snow, 920=Ice
            method = 1;
            
            k_vec = thermal_conductivity(t_vec, d_vec, tcase.rho_ice, method);
            
            % Expected Snow (Sturm)
            exp_snow = 0.138 - 1.01e-3 * d_vec(1) + 3.233e-6 * d_vec(1)^2;
            
            % Expected Ice (Temp dependent)
            exp_ice = 9.828 * exp(-5.7e-3 * t_vec(2));
            
            tcase.verifyEqual(k_vec(1), exp_snow, 'AbsTol', 1e-8);
            tcase.verifyEqual(k_vec(2), exp_ice, 'AbsTol', 1e-8);
        end
        
        function test_density_threshold_boundary(tcase)
            % Test the exact boundary of density_ice
            % Code uses: sfIdx = d < density_ice - 1e-11
            
            d_vals = [tcase.rho_ice - 1e-10; tcase.rho_ice];
            t_val = 260;
            
            k_out = thermal_conductivity([t_val; t_val], d_vals, tcase.rho_ice, 1);
            
            % 1. d = rho - 1e-10. This is practically rho, but mathematically < rho.
            % However, the tolerance is 1e-11.
            % (rho - 1e-10) < (rho - 1e-11) is TRUE. So it should be SNOW.
            exp_snow = 0.138 - 1.01e-3 * d_vals(1) + 3.233e-6 * d_vals(1)^2;
            tcase.verifyEqual(k_out(1), exp_snow, 'AbsTol', 1e-8, 'Just below threshold should be snow');
            
            % 2. d = rho. This is > (rho - 1e-11). So it should be ICE.
            exp_ice = 9.828 * exp(-5.7e-3 * t_val);
            tcase.verifyEqual(k_out(2), exp_ice, 'AbsTol', 1e-8, 'At threshold should be ice');
        end
    end
end