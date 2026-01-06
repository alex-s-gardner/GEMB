classdef test_layer_management < matlab.unittest.TestCase
    
    properties
        % Grid state vectors
        n = 10;
        t_vec
        dz
        d
        w
        re
        gdn
        gsp
        a
        a_diff
        
        % Configuration parameters
        dz_min = 0.05;      % Min layer thickness
        z_max_total = 100;  % Max column depth
        z_min_total = 0.5;  % Min column depth (Set LOW to avoid auto-adding layers)
        z_top = 2;          % Depth of constant grid spacing
        z_y = 1.1;          % Stretching factor
        verbose = false;
        
        % Constants from source
        ci = 2102;
        lf = 0.3345E6;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize a standard column
            % Total depth = 10 * 0.1 = 1.0m. 
            % This is > z_min_total (0.5), so no layers will be added automatically.
            tcase.t_vec = 260 * ones(tcase.n, 1);
            tcase.dz = 0.1 * ones(tcase.n, 1); 
            tcase.d = 400 * ones(tcase.n, 1);
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.5 * ones(tcase.n, 1);
            tcase.gdn = 0.5 * ones(tcase.n, 1);
            tcase.gsp = 0.5 * ones(tcase.n, 1);
            tcase.a = 0.8 * ones(tcase.n, 1);
            tcase.a_diff = 0.8 * ones(tcase.n, 1);
            
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
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, m_add, e_add] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n);
            tcase.verifyEqual(dz_out, tcase.dz);
            tcase.verifyEqual(t_out, tcase.t_vec);
            tcase.verifyEqual(m_add, 0);
            tcase.verifyEqual(e_add, 0);
        end
        
        function test_merge_small_layer(tcase)
            % Case: Top layer is too small (0.01 < 0.05). Should merge into layer 2.
            
            tcase.dz(1) = 0.01;
            
            % CRITICAL FIX: Reduce neighbor size to prevent immediate re-splitting.
            % dz_min = 0.05, so dz_max = 0.10.
            % If Layer 2 is 0.1, Sum = 0.11 -> Splits back to 0.055.
            % Set Layer 2 to 0.05. Sum = 0.06 (Valid).
            tcase.dz(2) = 0.05; 
            
            % Calculate expected mass weighted average for density/temp
            m1 = tcase.dz(1) * tcase.d(1);
            m2 = tcase.dz(2) * tcase.d(2);
            m_total = m1 + m2;
            expected_dz = tcase.dz(1) + tcase.dz(2);
            expected_d = m_total / expected_dz;
            
            [~, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~, ~] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n - 1, 'Should have 1 fewer layer');
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Layers should sum thickness');
            tcase.verifyEqual(d_out(1), expected_d, 'AbsTol', 1e-10, 'Density should be mass weighted');
        end
        
        function test_split_large_layer(tcase)
            % Case: Layer is too thick (> 2*dz_min or > stretching limit).
            % dz_min = 0.05. Max = 0.1 (in top section).
            % We set layer 1 to 0.2. It should split into two 0.1 layers.
            
            tcase.dz(1) = 0.2; 
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~, ~] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Should have 1 extra layer');
            tcase.verifyEqual(dz_out(1), 0.1, 'AbsTol', 1e-10, 'Layer 1 should be halved');
            tcase.verifyEqual(dz_out(2), 0.1, 'AbsTol', 1e-10, 'Layer 2 should be halved');
            tcase.verifyEqual(t_out(1), t_out(2), 'Split layers should share temp');
            tcase.verifyEqual(d_out(1), d_out(2), 'Split layers should share density');
        end
        
        function test_add_bottom_layer(tcase)
            % Case: Total column depth < z_min_total.
            % Code duplicates bottom layer to extend depth.
            
            % Setup a very shallow column
            tcase.dz = 0.1 * ones(5, 1); % Total 0.5m
            min_req = 1.0; % Requirement is 1.0m
            
            [t_out, dz_out, ~, ~, ~, ~, ~, ~, ~, m_add, ~] = layer_management(...
                tcase.t_vec(1:5), tcase.dz, tcase.d(1:5), tcase.w(1:5), tcase.re(1:5), ...
                tcase.gdn(1:5), tcase.gsp(1:5), tcase.a(1:5), tcase.a_diff(1:5), ...
                tcase.dz_min, tcase.z_max_total, min_req, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), 6, 'Should add 1 layer to meet depth req');
            tcase.verifyEqual(dz_out(end), dz_out(end-1), 'New layer should match previous bottom');
            
            % Check return variable for added mass
            m_expected = tcase.dz(end) * tcase.d(end);
            tcase.verifyEqual(m_add, m_expected, 'AbsTol', 1e-10, 'Should report added mass');
        end
        
        function test_remove_bottom_layer(tcase)
            % Case: Total column depth > z_max_total.
            % Code removes bottom layer.
            
            tcase.dz = 0.1 * ones(10, 1); 
            max_allowed = 0.95; 
            
            [~, dz_out, ~, ~, ~, ~, ~, ~, ~, m_add, ~] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, max_allowed, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), 9, 'Should remove bottom layer');
            
            % Check return variable for removed mass (negative)
            m_removed_expected = -(0.1 * 400); 
            tcase.verifyEqual(m_add, m_removed_expected, 'AbsTol', 1e-10, 'Should report removed mass (negative)');
        end
        
        function test_bottom_temp_boundary_condition(tcase)
            % The function enforces T(end) = T_bottom_initial to stabilize diffusion.
            % Verify this reset happens even if layers change.
            
            t_orig_bottom = tcase.t_vec(end);
            
            % Force a split at the top so indices shift, ensuring logic relies on values not just indices
            tcase.dz(1) = 0.2; 
            
            [t_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(t_out(end), t_orig_bottom, 'Bottom temperature must be clamped to original bottom temp');
        end
        
        function test_conservation_check(tcase)
            % Run with verbose=true to trigger the internal mass/energy check.
            % We create a scenario with complex merges and splits.
            
            tcase.verbose = true;
            tcase.dz(1) = 0.01; % Forces merge
            % Ensure neighbor is small enough to merge without splitting again
            tcase.dz(2) = 0.05; 
            
            try
                layer_management(...
                    tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                    tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                    tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                    tcase.z_top, tcase.z_y, tcase.verbose);
            catch ME
                tcase.verifyFail(['Conservation check failed: ' ME.message]);
            end
        end
        
        function test_bottom_merge_logic(tcase)
            % Special case: Bottom layer (index m) is too small.
            % Logic: "i_target = find(~delete_cell,1,'last');" -> Merges into the one above it.
            
            tcase.dz(end) = 0.01; % Bottom too small
            
            % CRITICAL FIX: Reduce neighbor (2nd to last) so combined size < 0.10.
            tcase.dz(end-1) = 0.05;
            
            [~, dz_out, ~, ~, ~, ~, ~, ~, ~, ~, ~] = layer_management(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.dz_min, tcase.z_max_total, tcase.z_min_total, ...
                tcase.z_top, tcase.z_y, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n - 1, 'Bottom layer should merge up');
            
            expected_dz_last = tcase.dz(end) + tcase.dz(end-1);
            tcase.verifyEqual(dz_out(end), expected_dz_last, 'AbsTol', 1e-10, 'Last layer should be sum of bottom two');
        end
    end
end