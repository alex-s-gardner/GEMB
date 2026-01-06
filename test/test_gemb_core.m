classdef test_gemb_core < matlab.unittest.TestCase
    
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
        
        % Forcing data
        dt = 3600; % 1 hour
        p = 0;
        ec = 0;
        m_surf = 0;
        density_ice = 917;
        bc_snow = 0;
        bc_ice = 0;
        sza = 60;
        cot = 0;
        cld_frac = 0;
        dsw = 200;
        dsw_diff = 50;
        dlw = 300;
        t_air = 265;
        v_air = 5;
        e_air = 400;
        p_air = 100000;
        
        % Settings Struct (S)
        S
        verbose = false;
    end
    
    methods (TestClassSetup)
        function create_mocks(tcase)
            % Add source path
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
        end
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize standard profile
            tcase.t_vec = 260 * ones(tcase.n, 1);
            
            % CRITICAL FIX: Set dz to 0.08.
            % dz_min = 0.05, so dz_max = 0.10.
            % Setting dz=0.10 puts it right on the splitting threshold.
            % 0.08 is safely in the middle, preventing unexpected splits.
            tcase.dz = 0.08 * ones(tcase.n, 1);
            
            tcase.d = 400 * ones(tcase.n, 1);
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.5 * ones(tcase.n, 1);
            tcase.gdn = 0.5 * ones(tcase.n, 1);
            tcase.gsp = 0.5 * ones(tcase.n, 1);
            tcase.a = 0.8 * ones(tcase.n, 1);
            tcase.a_diff = 0.8 * ones(tcase.n, 1);
            
            % Initialize Settings Struct S with reasonable defaults
            tcase.S.albedo_method = 1;
            tcase.S.albedo_ice = 0.45;
            tcase.S.albedo_snow = 0.85;
            tcase.S.albedo_fixed = 0.7;
            tcase.S.albedo_desnity_threshold = 1023;
            tcase.S.albedo_wet_snow_t0 = 15;
            tcase.S.albedo_dry_snow_t0 = 30;
            tcase.S.albedo_K = 7;
            tcase.S.sw_absorption_method = 1;
            tcase.S.Vz = 2;
            tcase.S.Tz = 2;
            tcase.S.emissivity = 0.98;
            tcase.S.ulw_delta = 0;
            tcase.S.emissivity_re_threshold = 10;
            tcase.S.emissivity_method = 0;
            tcase.S.thermal_conductivity_method = 1;
            tcase.S.T_mean = 260;
            tcase.S.column_dzmin = 0.05;
            tcase.S.P_mean = 200;
            tcase.S.V_mean = 5;
            tcase.S.new_snow_method = 1;
            tcase.S.column_zmax = 100;
            
            % Set zmin low to prevent auto-add. 
            % Current depth = 10 * 0.08 = 0.8m. zmin=0.5m is safe.
            tcase.S.column_zmin = 0.5; 
            
            tcase.S.column_ztop = 2;
            tcase.S.column_zy = 1.1;
            tcase.S.densification_method = 1;
        end
    end
    
    methods (Test)
        
        function test_pipeline_execution(tcase)
            % Basic "Smoke Test": Verify it runs without error and returns valid shapes
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, a_out, ~, ~, ~, sw_net, ...
             shf, lhf, ulw, ~, m_tot, r_tot, f_tot, m_add, e_add, comp_dens, comp_melt] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.dt, tcase.p, ...
                tcase.ec, tcase.m_surf, tcase.density_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, tcase.dsw, tcase.dsw_diff, ...
                tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, ...
                tcase.S, tcase.verbose);
            
            % Check Outputs
            tcase.verifyEqual(size(t_out), size(tcase.t_vec));
            tcase.verifyEqual(size(dz_out), size(tcase.dz));
            tcase.verifyEqual(size(d_out), size(tcase.d));
            
            % Physics Checks (Basic)
            tcase.verifyTrue(all(d_out > 0), 'Density must be positive');
            tcase.verifyTrue(all(t_out > 0), 'Temperature must be positive (Kelvin)');
            tcase.verifyTrue(all(a_out >= 0 & a_out <= 1), 'Albedo must be [0,1]');
            
            % Fluxes
            tcase.verifyTrue(isscalar(sw_net));
            tcase.verifyTrue(isscalar(shf));
            tcase.verifyTrue(isscalar(lhf));
            tcase.verifyTrue(isscalar(ulw));
            
            % Mass Balance Terms
            tcase.verifyTrue(isscalar(m_tot));
            tcase.verifyTrue(isscalar(r_tot));
            tcase.verifyTrue(isscalar(f_tot));
            tcase.verifyTrue(isscalar(m_add));
            tcase.verifyTrue(isscalar(e_add));
            tcase.verifyTrue(isscalar(comp_dens));
            tcase.verifyTrue(isscalar(comp_melt));
        end
        
        function test_accumulation_event(tcase)
            % Test that precipitation adds a layer (via accumulation sub-function)
            
            % Use 20 kg/m2 precip.
            % 20 kg/m2 / 350 kg/m3 = ~0.057m.
            % 0.057 > dz_min (0.05) -> NEW LAYER created.
            % 0.057 < 2*dz_min (0.10) -> NO SPLIT happens.
            tcase.p = 20; 
            tcase.t_air = 260; % Snow
            
            [~, dz_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.dt, tcase.p, ...
                tcase.ec, tcase.m_surf, tcase.density_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, tcase.dsw, tcase.dsw_diff, ...
                tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, ...
                tcase.S, tcase.verbose);
            
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Precipitation should add exactly 1 layer');
        end
        
        function test_melt_event(tcase)
            % Test that high temperature causes melt (via thermo -> melt sub-functions)
            tcase.t_air = 280; % Hot air
            tcase.dsw = 800;   % Intense sun
            tcase.t_vec(1) = 273.15; % Surface at melting point
            
            % Increase timestep to allow significant energy input
            tcase.dt = 3600 * 3; 
            
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, m_tot, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.dt, tcase.p, ...
                tcase.ec, tcase.m_surf, tcase.density_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, tcase.dsw, tcase.dsw_diff, ...
                tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, ...
                tcase.S, tcase.verbose);
            
            tcase.verifyTrue(m_tot > 0, 'High energy input should generate melt');
        end
        
        function test_densification_compaction(tcase)
            % Test that densification returns positive compaction value
            tcase.S.densification_method = 1; 
            
            % Set density low so there is room to compact
            tcase.d(:) = 300; 
            
            [~, ~, d_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, comp_dens, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.dt, tcase.p, ...
                tcase.ec, tcase.m_surf, tcase.density_ice, tcase.bc_snow, tcase.bc_ice, ...
                tcase.sza, tcase.cot, tcase.cld_frac, tcase.dsw, tcase.dsw_diff, ...
                tcase.dlw, tcase.t_air, tcase.v_air, tcase.e_air, tcase.p_air, ...
                tcase.S, tcase.verbose);
            
            tcase.verifyTrue(comp_dens > 0, 'Densification should result in positive dry compaction');
            tcase.verifyTrue(all(d_out >= 300), 'Density should increase');
        end
    end
end