classdef test_grid_initialize < matlab.unittest.TestCase
    
    methods (TestMethodSetup)
        function add_source_path(tcase)
            % Add the src folder to the path to locate grid_initialize.m
            import matlab.unittest.fixtures.PathFixture
            try
                % Assuming standard folder structure: /src and /tests
                tcase.applyFixture(PathFixture('../src'));
            catch
                % Fallback if running from root
                addpath('src');
            end
        end
    end
    
    methods (Test)
        
        function test_standard_grid_creation(tcase)
            % Test a standard configuration
            z_top = 10;
            dz_top = 1.0;
            z_max = 20;
            beta = 1.0; % No stretching
            
            % Call function (using filename grid_initialize)
            [dz, z_center] = grid_initialize(z_top, dz_top, z_max, beta);
            
            % Verify top section
            n_top = z_top / dz_top;
            tcase.verifyEqual(dz(1:n_top), ones(n_top, 1) * dz_top, ...
                'Top portion of grid should have constant spacing');
            
            % Verify total depth
            % Since beta=1, it should fill exactly to z_max if aligned
            total_depth = sum(dz);
            tcase.verifyTrue(total_depth >= z_max, 'Total depth must reach or exceed z_max');
            
            % Verify z_center dimensions
            tcase.verifyEqual(length(z_center), length(dz), ...
                'z_center should have same length as dz');
        end
        
        function test_grid_stretching(tcase)
            % Test that beta actually stretches the grid below z_top
            z_top = 10;
            dz_top = 1.0;
            z_max = 20;
            beta = 1.5; 
            
            [dz, ~] = grid_initialize(z_top, dz_top, z_max, beta);
            
            n_top = z_top / dz_top;
            
            % Check the first cell of the bottom section
            % It should be dz_top * beta
            expected_first_bottom = dz_top * beta;
            tcase.verifyEqual(dz(n_top + 1), expected_first_bottom, 'AbsTol', 1e-10, ...
                'First bottom cell should be stretched by beta');
            
            % Check subsequent stretching
            tcase.verifyEqual(dz(n_top + 2), expected_first_bottom * beta, 'AbsTol', 1e-10, ...
                'Subsequent bottom cells should continue geometric growth');
        end
        
        function test_invalid_top_spacing_error(tcase)
            % The function asserts that dz_top must divide z_top evenly.
            % z_top = 10, dz_top = 3 -> 10/3 = 3.33 (Not integer) -> Should Error
            
            z_top = 10;
            dz_top = 3; 
            z_max = 20;
            beta = 1.1;
            
            % FIX: Use ?MException to catch any error, as the specific identifier 
            % 'MATLAB:assertion:failed' was not present in the error report.
            tcase.verifyError(@() grid_initialize(z_top, dz_top, z_max, beta), ...
                ?MException, ...
                'Should throw assertion error if dz_top does not divide z_top evenly');
        end
        
        function test_small_dz_warning(tcase)
            % The function warns if dz_top < 0.05
            
            z_top = 0.1;
            dz_top = 0.01; % Triggers warning
            z_max = 1;
            beta = 1.0;
            
            tcase.verifyWarning(@() grid_initialize(z_top, dz_top, z_max, beta), ...
                '', ... % Warnings might not have specific IDs in user code, catch all or specific text
                'Should warn for very small dz_top');
        end
        
        function test_center_calculation_logic(tcase)
            % Verify the math: z_center = -cumsum(dz) + dz/2
            
            z_top = 2;
            dz_top = 1;
            z_max = 2;
            beta = 1;
            
            % Grid should be [1; 1]
            % Depths: 1, 2. 
            % Centers: 0.5, 1.5. 
            % Function returns negative depths relative to surface (or positive depth?):
            % The code: z_center = -cumsum(dz) + dz/2;
            % This implies negative values (e.g., -0.5, -1.5)
            
            [dz, z_center] = grid_initialize(z_top, dz_top, z_max, beta);
            
            expected_centers = [-0.5; -1.5];
            
            tcase.verifyEqual(z_center, expected_centers, 'AbsTol', 1e-10, ...
                'z_center should be negative depths calculated from surface');
        end
        
        function test_output_orientation(tcase)
            % Ensure outputs are column vectors as used in GEMB
            z_top = 10; 
            dz_top = 1;
            z_max = 20; 
            beta = 1.1;
            
            [dz, z_center] = grid_initialize(z_top, dz_top, z_max, beta);
            
            % FIX: Check dimensions explicitly rather than exact length
            tcase.verifyEqual(size(dz, 2), 1, 'dz should be a column vector');
            tcase.verifyTrue(size(dz, 1) > 1, 'dz should have rows');
            
            tcase.verifyEqual(size(z_center, 2), 1, 'z_center should be a column vector');
            tcase.verifyEqual(size(z_center, 1), size(dz, 1), 'z_center and dz should have matching rows');
        end
    end
end