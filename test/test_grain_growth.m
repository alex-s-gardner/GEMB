classdef test_grain_growth < matlab.unittest.TestCase
    
    properties
        % tolerance for floating point comparisons
        tol = 1e-6;
    end
    
    methods (TestMethodSetup)
        function add_source_path(tcase)
            % add the src folder to path so the function is found
            % assumes structure:
            %   /src/grain_growth.m
            %   /tests/test_grain_growth.m
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                % fallback if running from root
                addpath('src');
            end
        end
    end
    
    methods (Test)
        
        function test_albedo_method_skip(tcase)
            % verify that albedo_method values other than 1 or 2 return inputs unchanged
            n = 5;
            t_in = 260 * ones(n, 1);
            dz = 0.1 * ones(n, 1);
            d = 300 * ones(n, 1);
            w = zeros(n, 1);
            re_in = 0.5 * ones(n, 1);
            gdn_in = 0.5 * ones(n, 1);
            gsp_in = 0.5 * ones(n, 1);
            dt = 86400;
            
            % test method 0 (should skip)
            [re_out, gdn_out, gsp_out] = grain_growth(t_in, dz, d, w, re_in, gdn_in, gsp_in, dt, 0);
            tcase.verifyEqual(re_out, re_in);
            tcase.verifyEqual(gdn_out, gdn_in);
            tcase.verifyEqual(gsp_out, gsp_in);

            % test method 3 (should skip)
            [re_out, gdn_out, gsp_out] = grain_growth(t_in, dz, d, w, re_in, gdn_in, gsp_in, dt, 3);
            tcase.verifyEqual(re_out, re_in);
        end

        function test_dendritic_dry_low_gradient(tcase)
            % case: dendritic snow (gdn > 0), dry (w=0), low gradient (dt < 5)
            % setup 3 layers with uniform t to ensure gradient is 0
            
            method = 1;
            dt = 86400; % 1 day
            
            t_vec = [260; 260; 260]; 
            dz = [0.1; 0.1; 0.1];
            d = [200; 200; 200];
            w = [0; 0; 0];
            re_in = [0.2; 0.2; 0.2];
            gdn_in = [0.8; 0.8; 0.8]; % high dentricity
            gsp_in = [0.2; 0.2; 0.2];
            
            [re_out, gdn_out, gsp_out] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % expect gdn to decrease (decay) and gsp to increase
            tcase.verifyTrue(all(gdn_out < gdn_in), 'dentricity should decrease for dry low gradient');
            tcase.verifyTrue(all(gsp_out > gsp_in), 'sphericity should increase for dry low gradient');
            
            % re is updated based on gdn/gsp for dendritic snow
            tcase.verifyNotEqual(re_out, re_in);
        end

        function test_dendritic_dry_high_gradient(tcase)
            % case: dendritic snow, dry, high gradient (> 5 k/m)
            % setup t to create gradient. 
            % z centers approx: 0.05, 0.15, 0.25. dist = 0.1m.
            % delta t = 10 k. gradient = 100 k/m > 5.
            
            method = 1;
            dt = 86400;
            
            t_vec = [260; 270; 280]; 
            dz = [0.1; 0.1; 0.1];
            d = [200; 200; 200];
            w = [0; 0; 0];
            re_in = [0.2; 0.2; 0.2];
            gdn_in = [0.8; 0.8; 0.8];
            gsp_in = [0.2; 0.2; 0.2];
            
            [~, gdn_out, gsp_out] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % under high gradient, the coefficient c (negative) is applied to both
            tcase.verifyTrue(all(gdn_out < gdn_in), 'dentricity should decrease under high gradient');
            tcase.verifyTrue(all(gsp_out < gsp_in), 'sphericity should decrease under high gradient');
        end

        function test_dendritic_wet_snow(tcase)
            % case: dendritic snow, wet (w > 0)
            
            method = 1;
            dt = 86400;
            
            t_vec = [273.15; 273.15; 273.15];
            dz = [0.1; 0.1; 0.1];
            d = [250; 250; 250];
            % add water content
            w = [1.0; 1.0; 1.0]; 
            re_in = [0.2; 0.2; 0.2];
            gdn_in = [0.8; 0.8; 0.8];
            gsp_in = [0.2; 0.2; 0.2];
            
            [~, gdn_out, gsp_out] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % wet snow causes rapid rounding: gdn decreases, gsp increases
            tcase.verifyTrue(all(gdn_out < gdn_in), 'wet snow should reduce dentricity');
            tcase.verifyTrue(all(gsp_out > gsp_in), 'wet snow should increase sphericity');
        end

        function test_nondendritic_dry_marbouty(tcase)
            % case: gdn = 0, dry. marbouty model active.
            
            method = 1;
            dt = 86400;
            
            % gradient needed for marbouty g factor. use moderate gradient ~ 20 k/m
            t_vec = [250; 252; 254]; 
            dz = [0.1; 0.1; 0.1];
            % density must be < 400 for growth in marbouty
            d = [300; 300; 300]; 
            w = [0; 0; 0];
            re_in = [0.5; 0.5; 0.5];
            gdn_in = [0; 0; 0];
            gsp_in = [0.5; 0.5; 0.5];
            
            [re_out, gdn_out, ~] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % expect growth in grain size
            tcase.verifyTrue(all(re_out > re_in), 'grain size should increase (marbouty)');
            tcase.verifyEqual(gdn_out, gdn_in); % should stay 0
        end

        function test_marbouty_density_limit(tcase)
            % case: marbouty h factor becomes 0 if density > 400.
            
            method = 1;
            dt = 86400;
            
            t_vec = [250; 252; 254]; 
            dz = [0.1; 0.1; 0.1];
            d = [450; 450; 450]; % > 400 threshold
            w = [0; 0; 0];
            re_in = [0.5; 0.5; 0.5];
            gdn_in = [0; 0; 0];
            gsp_in = [0.5; 0.5; 0.5];
            
            [re_out, ~, ~] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % expect no growth
            tcase.verifyEqual(re_out, re_in, 'AbsTol', 1e-10, 'no growth expected above density threshold');
        end
        
        function test_nondendritic_wet_snow(tcase)
            % case: wet snow (w > 0), non-dendritic (gdn=0). brun 1989.
            
            method = 1;
            dt = 86400;
            
            t_vec = [273.15; 273.15; 273.15]; 
            dz = [0.1; 0.1; 0.1];
            d = [300; 300; 300];
            % create significant liquid water content. 
            w = [1.5; 1.5; 1.5];
            
            re_in = [0.5; 0.5; 0.5];
            gdn_in = [0; 0; 0];
            gsp_in = [0.5; 0.5; 0.5];
            
            [re_out, gdn_out, ~] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % expect growth via 'e' factor
            tcase.verifyTrue(all(re_out > re_in), 'wet snow should grow grains (brun)');
            tcase.verifyEqual(gdn_out, gdn_in);
        end

        function test_clamping_limits(tcase)
            % verify that dentricity and sphericity are clamped to [0,1]
            
            method = 1;
            dt = 86400 * 100; % long time to force overshoot
            
            t_vec = [260; 260; 260];
            dz = [0.1; 0.1; 0.1];
            d = [200; 200; 200];
            w = [0; 0; 0];
            re_in = [0.2; 0.2; 0.2];
            gdn_in = [0.1; 0.1; 0.1]; % close to 0
            gsp_in = [0.2; 0.2; 0.2];
            
            [~, gdn_out, gsp_out] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            tcase.verifyTrue(all(gdn_out >= 0), 'gdn should have lower bound of 0');
            tcase.verifyTrue(all(gsp_out <= 1), 'gsp should have upper bound of 1');
        end

        function test_grain_size_cap(tcase)
            % verify grain size cap at 2mm for fully spherical grains (gsp=1)
            
            method = 1;
            dt = 86400 * 50; % enough time to grow large
            
            t_vec = [273.15; 273.15; 273.15];
            dz = [0.1; 0.1; 0.1];
            d = [300; 300; 300];
            w = [2.0; 2.0; 2.0]; % wet snow grows fast
            
            % start near the cap
            re_in = [1.9; 1.9; 1.9]; % radius 1.9mm (grain size 3.8mm)
            gdn_in = [0; 0; 0];
            gsp_in = [1; 1; 1]; % sphericity 1 -> cap is 2mm size (1mm radius)
            
            [re_out, ~, ~] = grain_growth(t_vec, dz, d, w, re_in, gdn_in, gsp_in, dt, method);
            
            % the code caps gsz at 2 if gsp=1. re = gsz/2. so re should be capped at 1.0.
            tcase.verifyTrue(all(re_out <= 1.0 + 1e-10), 'radius should be capped at 1.0mm (size 2mm) for spherical grains');
        end
    end
end