classdef test_model_initialize_column < matlab.unittest.TestCase
    
    properties
        % Default ModelParam structure
        MP
    end
    
    methods (TestMethodSetup)
        function setup_environment(tcase)
            % Add the src folder to the path
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
            
            % Initialize default parameters
            tcase.MP.column_ztop  = 10;
            tcase.MP.column_dztop = 0.05;
            tcase.MP.column_zmax  = 20;
            tcase.MP.column_zy    = 1.0;
        end
    end
    
    methods (Test)
        
        function test_standard_grid_creation(tcase)
            % Test a standard configuration
            tcase.MP.column_ztop  = 10;
            tcase.MP.column_dztop = 0.05;
            tcase.MP.column_zmax  = 20;
            tcase.MP.column_zy    = 1.0; % No stretching
            
            ModelParam = model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop,...
            column_zmax=tcase.MP.column_zmax,...
            column_zy=tcase.MP.column_zy);
            ClimateForcing.temperature_air_mean = 253.15; % -20 C

            [~, dz] = model_initialize_column(ModelParam, ClimateForcing);
            z_center = dz2z(dz); 

            % Verify top section
            n_top = tcase.MP.column_ztop / tcase.MP.column_dztop;
            tcase.verifyEqual(dz(1:n_top), ones(n_top, 1) * tcase.MP.column_dztop, ...
                'Top portion of grid should have constant spacing');
            
            % Verify total depth
            total_depth = sum(dz);
            
            tolerance = 1e-10; 
            tcase.verifyTrue(total_depth >= (tcase.MP.column_zmax-tolerance), 'Total depth must reach or exceed z_max');

            % Verify z_center dimensions
            tcase.verifyEqual(length(z_center), length(dz), ...
                'z_center should have same length as dz');

        end
        
        function test_grid_stretching(tcase)
            % Test that column_zy actually stretches the grid below ztop
            tcase.MP.column_ztop  = 10;
            tcase.MP.column_dztop = 0.05;
            tcase.MP.column_zmax  = 20;
            tcase.MP.column_zy    = 1.5; 
            
            ModelParam = model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop,...
            column_zmax=tcase.MP.column_zmax,...
            column_zy=tcase.MP.column_zy);
            ClimateForcing.temperature_air_mean = 253.15; % -20 C

            [~, dz] = model_initialize_column(ModelParam, ClimateForcing);
            
            n_top = tcase.MP.column_ztop / tcase.MP.column_dztop;
            
            % Check the first cell of the bottom section
            % It should be dz_top * column_zy
            expected_first_bottom = tcase.MP.column_dztop * tcase.MP.column_zy;
            tcase.verifyEqual(dz(n_top + 1), expected_first_bottom, 'AbsTol', 1e-10, ...
                'First bottom cell should be stretched by column_zy');
            
            % Check subsequent stretching
            tcase.verifyEqual(dz(n_top + 2), expected_first_bottom * tcase.MP.column_zy, 'AbsTol', 1e-10, ...
                'Subsequent bottom cells should continue geometric growth');
        end
        
        function test_invalid_top_spacing_error(tcase)
            % The function asserts that dz_top must divide z_top evenly.
            tcase.MP.column_ztop = 10;
            tcase.MP.column_dztop = 3; % 10/3 is not integer
            
            tcase.verifyError(@() model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop), ...
                ?MException, ...
                'Should throw assertion error if dz_top does not divide z_top evenly');
        end
        
        function test_small_dz_warning(tcase)
            % The function warns if dz_top < 0.05
            tcase.MP.column_ztop  = 0.1;
            tcase.MP.column_dztop = 0.01; % Triggers warning
            tcase.MP.column_zmax  = 1;
            tcase.MP.column_zy    = 1.0;
            
            ClimateForcing.temperature_air_mean = 253.15; % -20 C
            ModelParam = model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop,...
            column_zmax=tcase.MP.column_zmax,...
            column_zy=tcase.MP.column_zy);
            tcase.verifyWarning(@() model_initialize_column(ModelParam, ClimateForcing), ...
                '', ... 
                'Should warn for very small dz_top');
        end
        
        function test_center_calculation_logic(tcase)
            % Verify the math: z_center = -cumsum(dz) + dz/2
            tcase.MP.column_ztop  = 2;
            tcase.MP.column_dztop = 0.05;
            tcase.MP.column_zmax  = 2;
            tcase.MP.column_zy    = 1;
            
            % Grid should be [1; 1]
            % Depths: 1, 2. 
            % Centers: 0.5, 1.5 (negative)
            
            ModelParam = model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop,...
            column_zmax=tcase.MP.column_zmax,...
            column_zy=tcase.MP.column_zy);
            ModelParam.column_dztop = 1; % override the dztop value bc 1 would not get past model_initialize_parameters error checks.  

            ClimateForcing.temperature_air_mean = 253.15; % -20 C

            [~, dz] = model_initialize_column(ModelParam, ClimateForcing);
            z_center = dz2z(dz); 

            
            expected_centers = [-0.5; -1.5];
            
            tcase.verifyEqual(z_center, expected_centers, 'AbsTol', 1e-10, ...
                'z_center should be negative depths calculated from surface');
        end
        
        function test_output_orientation(tcase)
            % Ensure outputs are column vectors as used in GEMB
            tcase.MP.column_ztop  = 10;
            tcase.MP.column_dztop = 0.05;
            tcase.MP.column_zmax  = 20;
            tcase.MP.column_zy    = 1.1;
            
            ModelParam = model_initialize_parameters(column_ztop=tcase.MP.column_ztop,...
            column_dztop=tcase.MP.column_dztop,...
            column_zmax=tcase.MP.column_zmax,...
            column_zy=tcase.MP.column_zy);
            ClimateForcing.temperature_air_mean = 253.15; % -20 C

            [~, dz] = model_initialize_column(ModelParam, ClimateForcing);
            z_center = dz2z(dz); 
            
            tcase.verifyEqual(size(dz, 2), 1, 'dz should be a column vector');
            tcase.verifyTrue(size(dz, 1) > 1, 'dz should have rows');
            
            tcase.verifyEqual(size(z_center, 2), 1, 'z_center should be a column vector');
            tcase.verifyEqual(size(z_center, 1), size(dz, 1), 'z_center and dz should have matching rows');
        end
    end
end