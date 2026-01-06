classdef test_albedo < matlab.unittest.TestCase
    
    properties
        % Common default inputs
        n = 10;
        t_vec
        dz
        d
        w
        re
        a_in
        a_diffuse_in
        dt = 10800; % 3 hours in seconds
        
        % Scalars
        p = 0;
        ec = 0;
        m_surf = 0;
        rho_ice = 917;
        bc_snow = 0;
        bc_ice = 0;
        sza = 60;
        cot = 0;
        cld_frac = 0;
        
        % Albedo params
        alb_ice = 0.45;
        alb_snow = 0.85;
        alb_fixed = 0.7;
        alb_dens_thresh = 1023; % standard threshold (density of water/ice limit)
        alb_wet_t0 = 15;
        alb_dry_t0 = 30;
        alb_k = 7;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize common profile vectors
            tcase.t_vec = 260 * ones(tcase.n, 1);
            tcase.dz = 0.1 * ones(tcase.n, 1);
            tcase.d = 350 * ones(tcase.n, 1); % Fresh snow density
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.2 * ones(tcase.n, 1); % 0.2 mm grain radius
            tcase.a_in = 0.8 * ones(tcase.n, 1);
            tcase.a_diffuse_in = 0.8 * ones(tcase.n, 1);
            
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
        
        function test_method_0_fixed_albedo(tcase)
            % Verify Method 0 returns the fixed value regardless of other inputs
            method = 0;
            
            [a_out, ~] = albedo(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            tcase.verifyEqual(a_out(1), tcase.alb_fixed, 'Method 0 should return fixed albedo');
        end
        
        function test_density_threshold_override(tcase)
            % If surface density > threshold, albedo should be fixed even if method is 1-4
            method = 1;
            
            % Set surface density to be very high (e.g., solid ice > threshold)
            % If threshold is 1023, we need higher? 
            % The code logic is: if (threshold - d) < tolerance => fixed.
            % So if d >= threshold.
            d_high = tcase.d;
            d_high(1) = 1025; 
            
            [a_out, ~] = albedo(tcase.t_vec, tcase.dz, d_high, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            tcase.verifyEqual(a_out(1), tcase.alb_fixed, 'Density above threshold should force fixed albedo');
        end
        
        function test_method_1_gardner(tcase)
            % Test Method 1 (Gardner & Sharp)
            method = 1;
            
            [a_out, a_diff_out] = albedo(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            % Check ranges [0, 1]
            tcase.verifyTrue(a_out(1) > 0 && a_out(1) <= 1, 'Method 1 albedo out of range');
            tcase.verifyTrue(a_diff_out(1) > 0 && a_diff_out(1) <= 1, 'Method 1 diffuse albedo out of range');
            
            % Gardner method usually produces high albedo for small grains
            tcase.verifyTrue(a_out(1) > 0.7, 'Small grain size should yield high albedo in Method 1');
        end
        
        function test_method_2_brun(tcase)
            % Test Method 2 (Brun et al.)
            method = 2;
            
            [a_out, ~] = albedo(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            tcase.verifyTrue(a_out(1) > 0 && a_out(1) <= 1, 'Method 2 albedo out of range');
        end
        
        function test_method_3_density_linear(tcase)
            % Test Method 3 (Linear interpolation by density)
            method = 3;
            
            % Case 1: Fresh snow density (300)
            d_fresh = tcase.d;
            d_fresh(1) = 300; 
            
            [a_fresh, ~] = albedo(tcase.t_vec, tcase.dz, d_fresh, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, 0.5, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            % Should be close to alb_snow (0.85)
            tcase.verifyEqual(a_fresh(1), tcase.alb_snow, 'AbsTol', 0.01, ...
                'Method 3 should return fresh snow albedo for density 300');
            
            % Case 2: Ice density (917)
            d_ice = tcase.d;
            d_ice(1) = tcase.rho_ice;
            
            [a_ice_out, ~] = albedo(tcase.t_vec, tcase.dz, d_ice, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, 0.5, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
             
            % Should be close to alb_ice (0.45)
            tcase.verifyEqual(a_ice_out(1), tcase.alb_ice, 'AbsTol', 0.01, ...
                'Method 3 should return ice albedo for ice density');
        end
        
        function test_method_4_decay(tcase)
            % Test Method 4 (Exponential Decay)
            method = 4;
            
            % Set previous albedo high, expect decay
            a_prev = tcase.a_in;
            a_prev(1) = 0.85;
            
            dt_days = 86400; % 1 day
            
            [a_out, ~] = albedo(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                a_prev, tcase.a_diffuse_in, dt_days, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            % Verify decay happened (new albedo < old albedo)
            tcase.verifyTrue(a_out(1) < a_prev(1), 'Method 4 should decay albedo over time');
            tcase.verifyTrue(a_out(1) >= tcase.alb_ice, 'Albedo should not decay below ice albedo');
            
            % Test fresh snow reset (Precipitation > 0)
            % Note: Method 4 logic "a(1) = albedo_snow - ..." happens if there is precipitation
            % Actually, in code: a(1) = albedo_snow - (albedo_snow - a(1)) * exp(-P./z_snow);
            p_heavy = 100; % Large precipitation event
            
            [a_reset, ~] = albedo(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                a_prev, tcase.a_diffuse_in, dt_days, p_heavy, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            tcase.verifyTrue(a_reset(1) > a_out(1), 'Heavy precipitation should reset/increase albedo');
        end
        
        function test_firn_transition_logic(tcase)
            % Test the logic block: "if (albedo_method < 3) ... (d(1) >= density_phc)"
            % This logic handles the transition between snow and ice albedos 
            % for densities between pore close-off (830) and ice (917).
            
            method = 1;
            d_firn = tcase.d;
            d_firn(1) = 875; % Between 830 and 917
            
            % Run albedo
            [a_out, ~] = albedo(tcase.t_vec, tcase.dz, d_firn, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                tcase.alb_dens_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
            
            % The code manually overrides the Gardner result with a linear interpolation
            % Check that result is within bounds [alb_ice, alb_snow]
            tcase.verifyTrue(a_out(1) >= 0.4, 'Firn albedo too low');
            tcase.verifyTrue(a_out(1) <= 0.8, 'Firn albedo too high');
        end
        
        function test_ice_surface_modification(tcase)
             % Test the logic where surface density is effectively ice (>= 910/917)
             % but still below threshold (if threshold is set high like 1023)
             
             method = 1;
             d_ice_surf = tcase.d;
             d_ice_surf(1) = 920; % > density_ice
             
             % Ensure threshold is higher so we enter the calculation block
             high_thresh = 1023;
             
             [a_out, ~] = albedo(tcase.t_vec, tcase.dz, d_ice_surf, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                tcase.m_surf, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                high_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
             
             % Should be calculated based on M_surf (melt)
             % With 0 melt, it should be close to albedo_ice_max (0.58 in code constant)
             tcase.verifyTrue(abs(a_out(1) - 0.58) < 0.05, ...
                 'Ice surface albedo with no melt should be near max ice albedo (0.58)');
             
             % With high melt
             m_high = 5000; % lots of accumulated melt
             [a_melt, ~] = albedo(tcase.t_vec, tcase.dz, d_ice_surf, tcase.w, tcase.re, ...
                tcase.a_in, tcase.a_diffuse_in, tcase.dt, tcase.p, tcase.ec, ...
                m_high, tcase.rho_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, ...
                method, tcase.alb_ice, tcase.alb_snow, tcase.alb_fixed, ...
                high_thresh, tcase.alb_wet_t0, tcase.alb_dry_t0, tcase.alb_k);
             
             % Should decay towards min ice albedo
             tcase.verifyTrue(a_melt(1) < a_out(1), 'Melt should reduce ice albedo');
        end
    end
end