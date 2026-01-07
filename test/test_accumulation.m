classdef test_accumulation < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 5;
        t_vec
        dz
        d
        w
        re
        gdn
        gsp
        a_in
        a_diff_in
        
        % Constants and Parameters
        rho_ice = 917;
        t_mean = 260;
        p_mean = 200;
        v_mean = 5;
        alb_snow = 0.85;
        alb_method = 1;
        dz_min = 0.05; % Minimum thickness to create new layer
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize standard profiles (cold firn column)
            tcase.t_vec = 260 * ones(tcase.n, 1);
            tcase.dz = 0.1 * ones(tcase.n, 1);
            tcase.d = 400 * ones(tcase.n, 1); 
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.5 * ones(tcase.n, 1);
            tcase.gdn = 0.5 * ones(tcase.n, 1);
            tcase.gsp = 0.5 * ones(tcase.n, 1);
            tcase.a_in = 0.7 * ones(tcase.n, 1);
            tcase.a_diff_in = 0.7 * ones(tcase.n, 1);
            
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
        
        function test_no_precipitation(tcase)
            % Verify that zero precipitation results in no changes
            p = 0;
            t_air = 270;
            v = 5;
            method = 0;
            
            [t_out, dz_out, d_out, w_out, ~, ~, ~, ~, ~, ra_out] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, method);
            
            tcase.verifyEqual(dz_out, tcase.dz);
            tcase.verifyEqual(d_out, tcase.d);
            tcase.verifyEqual(t_out, tcase.t_vec);
            tcase.verifyEqual(ra_out, 0);
        end
        
        function test_large_snow_event_new_layer(tcase)
            % Large snowfall (> dz_min) should create a new layer on top
            p = 50; % kg m-2 (roughly 15-30cm of snow depending on density)
            t_air = 260; % Cold air -> Snow
            v = 5;
            method = 0; % Default density = 150
            rho_snow = 150;
            
            [t_out, dz_out, d_out, ~, ~, gdn_out, gsp_out, a_out, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, method);
            
            % Check vector growth
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Large snow should add a layer');
            
            % Check top layer properties
            tcase.verifyEqual(d_out(1), rho_snow, 'Top layer should have fresh snow density');
            expected_dz = p / rho_snow;
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10);
            tcase.verifyEqual(t_out(1), t_air);
            tcase.verifyEqual(a_out(1), tcase.alb_snow);
            
            % Check default microstructure for new snow
            tcase.verifyEqual(gdn_out(1), 1.0); % Dendricity default
            tcase.verifyEqual(gsp_out(1), 0.5); % Sphericity default
        end
        
        function test_small_snow_event_merge(tcase)
            % Small snowfall (< dz_min) should merge with top layer
            % dz_min is 0.05. Snow density 150.
            % Limit is 0.05 * 150 = 7.5 kg.
            p = 2; % 2 kg m-2
            t_air = 260;
            v = 5;
            method = 0;
            rho_snow = 150;
            
            old_mass = tcase.d(1) * tcase.dz(1);
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, a_out, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, method);
            
            % Verify no new layer
            tcase.verifyEqual(length(dz_out), tcase.n, 'Small snow should merge');
            
            % Verify mass conservation and mixing
            new_mass = old_mass + p;
            expected_dz = tcase.dz(1) + p/rho_snow;
            expected_d = new_mass / expected_dz;
            
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Depth should increase');
            tcase.verifyEqual(d_out(1), expected_d, 'AbsTol', 1e-10, 'Density should decrease (mix with light snow)');
            
            % Verify Temperature weighting
            expected_t = (t_air * p + tcase.t_vec(1) * old_mass) / new_mass;
            tcase.verifyEqual(t_out(1), expected_t, 'AbsTol', 1e-10, 'Temperature should be weighted average');
            
            % Verify Albedo update
            expected_a = (tcase.alb_snow * p + tcase.a_in(1) * old_mass) / new_mass;
            tcase.verifyEqual(a_out(1), expected_a, 'AbsTol', 1e-10, 'Albedo should be weighted average');
        end
        
        function test_rain_event(tcase)
            % Rain (T_air > 0C) should add mass and heat to top layer, not thickness
            p = 10;
            t_air = 275; % > 273.15 -> Rain
            v = 0;
            method = 0;
            
            % Constants from function
            lf = 0.3345e6;
            ci = 2102;
            
            old_mass = tcase.d(1) * tcase.dz(1);
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, ra_out] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, method);
            
            % Rain Output flag
            tcase.verifyEqual(ra_out, p);
            
            % Mass update (thickness stays same, density increases for rain in this model)
            % The code: dz(1) stays same, d(1) = mass/dz(1)
            new_mass = old_mass + p;
            tcase.verifyEqual(d_out(1), new_mass/tcase.dz(1), 'AbsTol', 1e-10);
            tcase.verifyEqual(dz_out(1), tcase.dz(1), 'Thickness unchanged for rain unless ice density hit');
            
            % Temperature update (includes Latent Heat logic)
            % T(1) = (P *(T_air + LF/CI) + T(1) * mInit(1)) / mass;
            term_rain = p * (t_air + lf/ci);
            term_snow = tcase.t_vec(1) * old_mass;
            expected_t = (term_rain + term_snow) / new_mass;
            
            tcase.verifyEqual(t_out(1), expected_t, 'AbsTol', 1e-8, 'Temperature should account for rain latent heat');
        end
        
        function test_rain_density_cap(tcase)
            % If rain causes density to exceed ice density, it should clamp density and increase thickness
            p = 500; % Huge rain event
            t_air = 275;
            
            % Start with density near ice
            tcase.d(1) = 900; 
            tcase.dz(1) = 0.1; % Mass = 90
            
            % New mass = 590. If dz kept 0.1, density = 5900 (impossible).
            % Expect density -> rho_ice, dz -> increases.
            
            [~, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, 0, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, 0);
            
            tcase.verifyEqual(d_out(1), tcase.rho_ice, 'AbsTol', 1e-10, 'Density should be capped at ice density');
            
            total_mass = 90 + 500;
            expected_dz = total_mass / tcase.rho_ice;
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Thickness should expand to conserve mass at ice density');
        end
        
        function test_density_methods(tcase)
            % Test the different switches for new snow density
            p = 50; % Ensure new layer
            t_air = 250;
            v = 5;
            
            % Method 1: Antarctica (350)
            [~, ~, d1, ~, ~, ~, ~, ~, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, 1);
            tcase.verifyEqual(d1(1), 350, 'Method 1 should be 350');
            
            % Method 2: Greenland (315)
            [~, ~, d2, ~, ~, ~, ~, ~, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, 2);
            tcase.verifyEqual(d2(1), 315, 'Method 2 should be 315');
            
            % Method 3: Kaspers (Formula)
            % dSnow=(7.36e-2 + 1.06e-3*min(T_air_mean,273.15) + 6.69e-2*P_mean/1000. + 4.77e-3*V_mean)*1000.
            expected_3 = (7.36e-2 + 1.06e-3*tcase.t_mean + 6.69e-2*tcase.p_mean/1000 + 4.77e-3*tcase.v_mean)*1000;
            [~, ~, d3, ~, ~, ~, ~, ~, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, 3);
            tcase.verifyEqual(d3(1), expected_3, 'AbsTol', 1e-5, 'Method 3 calculation incorrect');
            
            % Method 4: Kuipers Munneke
            % dSnow = 481.0 + 4.834*(T_air_mean-273.15);
            expected_4 = 481.0 + 4.834*(tcase.t_mean - 273.15);
            [~, ~, d4, ~, ~, ~, ~, ~, ~, ~] = accumulation(...
                tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a_in, tcase.a_diff_in, ...
                t_air, p, v, tcase.rho_ice, tcase.t_mean, tcase.dz_min, ...
                tcase.p_mean, tcase.v_mean, tcase.alb_snow, tcase.alb_method, 4);
            tcase.verifyEqual(d4(1), expected_4, 'AbsTol', 1e-5, 'Method 4 calculation incorrect');
        end
    end
end