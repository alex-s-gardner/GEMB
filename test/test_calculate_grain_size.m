classdef test_calculate_grain_size < matlab.unittest.TestCase
    
    properties
        % tolerance for floating point comparisons
        tol = 1e-6;
        
        % Structures
        CF
        MP
    end
    
    methods (TestMethodSetup)
        function setup_structs(tcase)
            % Initialize default structures
            tcase.CF.dt = 86400; % 1 day in seconds
            tcase.MP.albedo_method = "GardnerSharp";
            
            % Add the src folder to path
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
        end
    end
    
    methods (Test)
        
        function test_albedo_method_skip(tcase)
            % verify that albedo_method values other than GardnerSharp/BrunLefebre return inputs unchanged
            n                   = 5;
            t_in                = 260 * ones(n, 1);
            dz                  = 0.1 * ones(n, 1);
            density             = 300 * ones(n, 1);
            water               = zeros(n, 1);
            grain_size_in       = 0.5 * ones(n, 1);
            grain_dendricity_in = 0.5 * ones(n, 1);
            grain_sphericity_in = 0.5 * ones(n, 1);
            
            % test method "None" (should skip)
            tcase.MP.albedo_method = "None";
            [grain_size_out, grain_dendricity_out, grain_sphericity_out] = calculate_grain_size(t_in, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            tcase.verifyEqual(grain_size_out, grain_size_in);
            tcase.verifyEqual(grain_dendricity_out, grain_dendricity_in);
            tcase.verifyEqual(grain_sphericity_out, grain_sphericity_in);

            % test method "GreuellKonzelmann" (should skip)
            tcase.MP.albedo_method = "GreuellKonzelmann";
            [grain_size_out, ~, ~] = calculate_grain_size(t_in, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            tcase.verifyEqual(grain_size_out, grain_size_in);
        end

        function test_dendritic_dry_low_gradient(tcase)
            % case: dendritic snow (gdn > 0), dry (water=0), low gradient (dt < 5)
            
            tcase.MP.albedo_method = "GardnerSharp";
            tcase.CF.dt = 86400;
            
            t_vec               = [260; 260; 260]; 
            dz                  = [0.1; 0.1; 0.1];
            density             = [200; 200; 200];
            water               = [  0;   0;   0];
            grain_size_in       = [0.2; 0.2; 0.2];
            grain_dendricity_in = [0.8; 0.8; 0.8]; % high dentricity
            grain_sphericity_in = [0.2; 0.2; 0.2];
            
            [grain_size_out, grain_dendricity_out, grain_sphericity_out] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % expect gdn to decrease (decay) and gsp to increase
            tcase.verifyTrue(all(grain_dendricity_out < grain_dendricity_in), 'dentricity should decrease for dry low gradient');
            tcase.verifyTrue(all(grain_sphericity_out > grain_sphericity_in), 'sphericity should increase for dry low gradient');
            
            % grain_radius is updated based on gdn/gsp for dendritic snow
            tcase.verifyNotEqual(grain_size_out, grain_size_in);
        end

        function test_dendritic_dry_high_gradient(tcase)
            % case: dendritic snow, dry, high gradient (> 5 k/m)
            
            tcase.MP.albedo_method = "GardnerSharp";
            
            t_vec               = [260; 270; 280]; 
            dz                  = [0.1; 0.1; 0.1];
            density             = [200; 200; 200];
            water               = [  0;   0;   0];
            grain_size_in       = [0.2; 0.2; 0.2];
            grain_dendricity_in = [0.8; 0.8; 0.8];
            grain_sphericity_in = [0.2; 0.2; 0.2];
            
            [~, grain_dendricity_out, grain_sphericity_out] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % under high gradient, the coefficient c (negative) is applied to both
            tcase.verifyTrue(all(grain_dendricity_out < grain_dendricity_in), 'dentricity should decrease under high gradient');
            tcase.verifyTrue(all(grain_sphericity_out < grain_sphericity_in), 'sphericity should decrease under high gradient');
        end

        function test_dendritic_wet_snow(tcase)
            % case: dendritic snow, wet (water > 0)
            
            tcase.MP.albedo_method = "GardnerSharp";
            
            t_vec               = [273.15; 273.15; 273.15];
            dz                  = [0.1; 0.1; 0.1];
            density             = [250; 250; 250];
            water               = [1.0; 1.0; 1.0]; 
            grain_size_in       = [0.2; 0.2; 0.2];
            grain_dendricity_in = [0.8; 0.8; 0.8];
            grain_sphericity_in = [0.2; 0.2; 0.2];
            
            [~, grain_dendricity_out, grain_sphericity_out] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % wet snow causes rapid rounding: gdn decreases, gsp increases
            tcase.verifyTrue(all(grain_dendricity_out < grain_dendricity_in), 'wet snow should reduce dentricity');
            tcase.verifyTrue(all(grain_sphericity_out > grain_sphericity_in), 'wet snow should increase sphericity');
        end

        function test_nondendritic_dry_marbouty(tcase)
            % case: gdn = 0, dry. marbouty model active.
            
            tcase.MP.albedo_method = "GardnerSharp";
            
            % gradient needed for marbouty g factor. use moderate gradient ~ 20 k/m
            t_vec               = [250; 252; 254]; 
            dz                  = [0.1; 0.1; 0.1];
            % density must be < 400 for growth in marbouty
            density             = [300; 300; 300]; 
            water               = [  0;   0;   0];
            grain_size_in       = [0.5; 0.5; 0.5];
            grain_dendricity_in = [  0;   0;   0];
            grain_sphericity_in = [0.5; 0.5; 0.5];
            
            [grain_size_out, grain_dendricity_out, ~] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % expect growth in grain size
            tcase.verifyTrue(all(grain_size_out > grain_size_in), 'grain size should increase (marbouty)');
            tcase.verifyEqual(grain_dendricity_out, grain_dendricity_in); % should stay 0
        end

        function test_marbouty_density_limit(tcase)
            % case: marbouty h factor becomes 0 if density > 400.
            
            tcase.MP.albedo_method = "GardnerSharp";
            
            t_vec               = [250; 252; 254]; 
            dz                  = [0.1; 0.1; 0.1];
            density             = [450; 450; 450]; % > 400 threshold
            water               = [  0;   0;   0];
            grain_size_in       = [0.5; 0.5; 0.5];
            grain_dendricity_in = [  0;   0;   0];
            grain_sphericity_in = [0.5; 0.5; 0.5];
            
            [grain_size_out, ~, ~] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % expect no growth
            tcase.verifyEqual(grain_size_out, grain_size_in, 'AbsTol', 1e-10, 'no growth expected above density threshold');
        end
        
        function test_nondendritic_wet_snow(tcase)
            % case: wet snow (water > 0), non-dendritic (gdn=0). brun 1989.
            
            tcase.MP.albedo_method = "GardnerSharp";
            
            t_vec               = [273.15; 273.15; 273.15]; 
            dz                  = [0.1; 0.1; 0.1];
            density             = [300; 300; 300];
            water               = [1.5; 1.5; 1.5]; 
            grain_size_in       = [0.5; 0.5; 0.5];
            grain_dendricity_in = [  0;   0;   0];
            grain_sphericity_in = [0.5; 0.5; 0.5];
            
            [grain_size_out, grain_dendricity_out, ~] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            % expect growth via 'e' factor
            tcase.verifyTrue(all(grain_size_out > grain_size_in), 'wet snow should grow grains (brun)');
            tcase.verifyEqual(grain_dendricity_out, grain_dendricity_in);
        end

        function test_clamping_limits(tcase)
            % verify that dentricity and sphericity are clamped to [0,1]
            
            tcase.MP.albedo_method = "GardnerSharp";
            tcase.CF.dt = 86400 * 100;
            
            t_vec               = [260; 260; 260];
            dz                  = [0.1; 0.1; 0.1];
            density             = [200; 200; 200];
            water               = [  0;   0;   0];
            grain_size_in       = [0.2; 0.2; 0.2];
            grain_dendricity_in = [0.1; 0.1; 0.1]; % close to 0
            grain_sphericity_in = [0.2; 0.2; 0.2];
            
            [~, grain_dendricity_out, grain_sphericity_out] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(grain_dendricity_out >= 0), 'gdn should have lower bound of 0');
            tcase.verifyTrue(all(grain_sphericity_out <= 1), 'gsp should have upper bound of 1');
        end

        function test_grain_size_cap(tcase)
            % verify grain size cap at 2mm for fully spherical grains (gsp=1)
            
            tcase.MP.albedo_method = "GardnerSharp";
            tcase.CF.dt = 86400 * 50;
            
            t_vec               = [273.15; 273.15; 273.15];
            dz                  = [0.1; 0.1; 0.1];
            density             = [300; 300; 300];
            water               = [2.0; 2.0; 2.0];
            grain_size_in       = [1.9; 1.9; 1.9]; % radius 1.9mm
            grain_dendricity_in = [  0;   0;   0];
            grain_sphericity_in = [  1;   1;   1];
            
            [grain_size_out, ~, ~] = calculate_grain_size(t_vec, dz, density, water, grain_size_in, grain_dendricity_in, grain_sphericity_in, tcase.CF, tcase.MP);
            
            tcase.verifyTrue(all(grain_size_out <= 1.0 + 1e-10), 'radius should be capped at 1.0mm (size 2mm) for spherical grains');
        end
    end
end