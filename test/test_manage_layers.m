classdef test_manage_layers < matlab.unittest.TestCase
    
    properties
        % Grid state vectors
        n = 10;
        t_vec
        dz
        density
        water
        grain_radius
        grain_dendricity
        grain_sphericity
        albedo
        albedo_diffuse
        
        % Structure
        MP % ModelParam
        
        verbose = false;
        
        % Constants from source
        ci = 2102;
        lf = 0.3345E6;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize a standard column
            % Total depth = 10 * 0.1 = 1.0m. 
            tcase.t_vec            = 260 * ones(tcase.n, 1);
            tcase.dz               = 0.1 * ones(tcase.n, 1); 
            tcase.density          = 400 * ones(tcase.n, 1);
            tcase.water            = zeros(tcase.n, 1);
            tcase.grain_radius     = 0.5 * ones(tcase.n, 1);
            tcase.grain_dendricity = 0.5 * ones(tcase.n, 1);
            tcase.grain_sphericity = 0.5 * ones(tcase.n, 1);
            tcase.albedo           = 0.8 * ones(tcase.n, 1);
            tcase.albedo_diffuse   = 0.8 * ones(tcase.n, 1);
            
            % Initialize ModelParam (MP)
            tcase.MP.column_dzmin  = 0.05;    % Min layer thickness
            tcase.MP.column_dzmax  = 0.10;    % Max layer thickness
            
            % CRITICAL FIX: Set zmax to exactly the initial depth (1.0m)
            % This prevents the code from automatically adding padding layers
            % unless the test explicitly changes this value.
            tcase.MP.column_zmax   = 1.0;   
            tcase.MP.column_zmin   = 0.5;     % Min column depth
            tcase.MP.column_ztop   = 2;       % Depth of constant grid spacing
            tcase.MP.column_zy     = 1.1;     % Stretching factor
            
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
        
        function test_no_action_needed(tcase)
            % Case: All layers are within valid bounds. No merge, no split.
            % zmax is 1.0, total depth is 1.0. Should remain 10 layers.
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, m_add, e_add] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n);
            tcase.verifyEqual(dz_out, tcase.dz);
            tcase.verifyEqual(t_out, tcase.t_vec);
            tcase.verifyEqual(m_add, 0);
            tcase.verifyEqual(e_add, 0);
        end
        
        function test_merge_small_layer(tcase)
            % Case: Top layer is too small. Should merge into layer 2.
            
            tcase.dz(1) = 0.01;
            tcase.dz(2) = 0.05; 
            
            % Total depth now: 0.01 + 0.05 + 0.8 = 0.86m
            % We must update zmax so the code doesn't try to pad it back to 1.0m
            tcase.MP.column_zmax = sum(tcase.dz);
            
            m1 = tcase.dz(1) * tcase.density(1);
            m2 = tcase.dz(2) * tcase.density(2);
            m_total = m1 + m2;
            expected_dz = tcase.dz(1) + tcase.dz(2);
            expected_d = m_total / expected_dz;
            
            [~, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~, ~] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n - 1, 'Should have 1 fewer layer');
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Layers should sum thickness');
            tcase.verifyEqual(d_out(1), expected_d, 'AbsTol', 1e-10, 'Density should be mass weighted');
        end
        
        function test_split_large_layer(tcase)
            % Case: Layer is too thick (> dz_max).
            
            tcase.dz(1) = 0.2; 
            % Total depth is now 1.1m (0.2 + 0.9). 
            % Update zmax so code doesn't prune the bottom layer.
            tcase.MP.column_zmax = 1.1;
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~, ~] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Should have 1 extra layer');
            tcase.verifyEqual(dz_out(1), 0.1, 'AbsTol', 1e-10, 'Layer 1 should be halved');
            tcase.verifyEqual(dz_out(2), 0.1, 'AbsTol', 1e-10, 'Layer 2 should be halved');
            tcase.verifyEqual(t_out(1), t_out(2), 'Split layers should share temp');
            tcase.verifyEqual(d_out(1), d_out(2), 'Split layers should share density');
        end
        
        function test_add_bottom_layer(tcase)
            % Case: Total column depth < zmax.
            
            tcase.dz = 0.1 * ones(5, 1); % Total 0.5m
            tcase.MP.column_zmax = 1.0;  % Target is 1.0m -> Adds padding
            
            [t_out, dz_out, ~, ~, ~, ~, ~, ~, ~, m_add, ~] = manage_layers(...
                tcase.t_vec(1:5), tcase.dz, tcase.density(1:5), tcase.water(1:5), tcase.grain_radius(1:5), ...
                tcase.grain_dendricity(1:5), tcase.grain_sphericity(1:5), tcase.albedo(1:5), tcase.albedo_diffuse(1:5), ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), 6, 'Should add 1 layer to meet depth req');
            tcase.verifyEqual(dz_out(end), dz_out(end-1), 'New layer should match previous bottom');
            
            m_expected = tcase.dz(end) * tcase.density(end);
            tcase.verifyEqual(m_add, m_expected, 'AbsTol', 1e-10, 'Should report added mass');
        end
        
        function test_remove_bottom_layer(tcase)
            % Case: Total column depth > column_zmax.
            
            tcase.dz = 0.1 * ones(10, 1); % Total 1.0m
            tcase.MP.column_zmax = 0.95;  % Max allowed is 0.95m -> Removes bottom
            
            [~, dz_out, ~, ~, ~, ~, ~, ~, ~, m_add, ~] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), 9, 'Should remove bottom layer');
            
            m_removed_expected = -(0.1 * 400); 
            tcase.verifyEqual(m_add, m_removed_expected, 'AbsTol', 1e-10, 'Should report removed mass (negative)');
        end
        
        function test_bottom_temp_boundary_condition(tcase)
            % Verify T(end) reset.
            
            t_orig_bottom = tcase.t_vec(end);
            tcase.dz(1) = 0.2; 
            tcase.MP.column_zmax = 1.1; % Adjust zmax to match new depth
            
            [t_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(t_out(end), t_orig_bottom, 'Bottom temperature must be clamped to original bottom temp');
        end
        
        function test_conservation_check(tcase)
            % Run with verbose=true to trigger the internal mass/energy check.
            
            tcase.verbose = true;
            tcase.dz(1) = 0.01; 
            tcase.dz(2) = 0.05; 
            tcase.MP.column_zmax = sum(tcase.dz); % Match depth
            
            try
                manage_layers(...
                    tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                    tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                    tcase.MP, tcase.verbose);
            catch ME
                tcase.verifyFail(['Conservation check failed: ' ME.message]);
            end
        end
        
        function test_bottom_merge_logic(tcase)
            % Bottom layer too small -> Merges up.
            
            tcase.dz(end) = 0.01; 
            tcase.dz(end-1) = 0.05;
            
            % Update zmax to match new total depth so padding doesn't happen
            tcase.MP.column_zmax = sum(tcase.dz);
            
            [~, dz_out, ~, ~, ~, ~, ~, ~, ~, ~, ~] = manage_layers(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, ...
                tcase.MP, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n - 1, 'Bottom layer should merge up');
            
            expected_dz_last = tcase.dz(end) + tcase.dz(end-1);
            tcase.verifyEqual(dz_out(end), expected_dz_last, 'AbsTol', 1e-10, 'Last layer should be sum of bottom two');
        end
    end
end