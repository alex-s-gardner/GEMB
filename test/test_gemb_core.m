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
        
        % State scalars
        EC_prev = 0;
        M_surf_prev = 0;
        
        % Structs
        CF % ClimateForcingStep
        MP % ModelParam
        
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
            
            % Initial depth: 10 layers * 0.08m = 0.8m Total Depth
            tcase.dz = 0.08 * ones(tcase.n, 1);
            
            tcase.d = 400 * ones(tcase.n, 1);
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.5 * ones(tcase.n, 1);
            tcase.gdn = 0.5 * ones(tcase.n, 1);
            tcase.gsp = 0.5 * ones(tcase.n, 1);
            tcase.a = 0.8 * ones(tcase.n, 1);
            tcase.a_diff = 0.8 * ones(tcase.n, 1);
            
            % --- Initialize ClimateForcingStep (CF) ---
            tcase.CF.dt = 3600; % 1 hour
            tcase.CF.P = 0;
            tcase.CF.T_air = 265;
            tcase.CF.V = 5;
            tcase.CF.e_air = 400;
            tcase.CF.p_air = 100000;
            tcase.CF.dlw = 300;
            tcase.CF.dsw = 200;
            tcase.CF.dsw_diff = 50;
            tcase.CF.black_carbon_snow = 0;
            tcase.CF.black_carbon_ice = 0;
            tcase.CF.solar_zenith_angle = 60;
            tcase.CF.cloud_optical_thickness = 0;
            tcase.CF.cloud_fraction = 0;
            
            % Location/Mean parameters
            tcase.CF.Vz = 2;
            tcase.CF.Tz = 2;
            tcase.CF.T_air_mean = 260;
            tcase.CF.P_mean = 200;
            tcase.CF.V_mean = 5;
            
            % --- Initialize ModelParam (MP) ---
            tcase.MP.albedo_method = "GardnerSharp"; 
            tcase.MP.albedo_ice = 0.45;
            tcase.MP.albedo_snow = 0.85;
            tcase.MP.albedo_fixed = 0.7;
            tcase.MP.albedo_desnity_threshold = 1023;
            tcase.MP.albedo_wet_snow_t0 = 15;
            tcase.MP.albedo_dry_snow_t0 = 30;
            tcase.MP.albedo_K = 7;
            tcase.MP.sw_absorption_method = 1;
            tcase.MP.emissivity = 0.98;
            tcase.MP.ulw_delta = 0;
            tcase.MP.emissivity_re_threshold = 10;
            tcase.MP.emissivity_method = 0;
            tcase.MP.thermal_conductivity_method = "Sturm"; 
            tcase.MP.column_dzmin = 0.05;
            tcase.MP.column_dzmax = 0.10; 
            tcase.MP.new_snow_method = "150kgm2"; 
            tcase.MP.density_ice = 917;
            
            % Set zmax slightly larger than initial depth to prevent removal.
            % Initial depth = 0.8m. 
            tcase.MP.column_zmax = 0.9; 
            tcase.MP.column_zmin = 0.5; 
            
            tcase.MP.column_ztop = 2;
            tcase.MP.column_zy = 1.1;
            tcase.MP.densification_method = "HerronLangway"; 
        end
    end
    
    methods (Test)
        
        function test_pipeline_execution(tcase)
            % Basic "Smoke Test": Verify it runs without error
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, a_out, ~, ~, ~, sw_net, ...
             shf, lhf, ulw, ~, m_tot, r_tot, f_tot, m_add, e_add, comp_dens, comp_melt] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.EC_prev, tcase.M_surf_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            % Check Outputs
            % Size logic: 
            % Initial 10 layers (0.8m). zmax=0.9m.
            % layer_management sees depth < zmax -> adds 1 layer.
            % Result: 11 layers.
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Expect 1 padding layer');
            
            tcase.verifyEqual(size(dz_out), size(t_out));
            tcase.verifyEqual(size(dz_out), size(d_out));
            
            % Fluxes
            tcase.verifyTrue(isscalar(sw_net));
            tcase.verifyTrue(isscalar(shf));
            tcase.verifyTrue(isscalar(lhf));
            tcase.verifyTrue(isscalar(ulw));
        end
        
        function test_accumulation_event(tcase)
            % Test that precipitation adds a layer
            
            % 1. Use P = 10 kg/m2.
            tcase.CF.P = 10; 
            tcase.CF.T_air = 260; % Snow
            
            % 2. Set zmax huge to force padding logic instead of removal logic.
            tcase.MP.column_zmax = 100.0;
            
            [~, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.EC_prev, tcase.M_surf_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            % Expected Layers:
            % 10 (Original) + 1 (Snow Accumulation) + 1 (Zmax Padding) = 12
            tcase.verifyEqual(length(dz_out), tcase.n + 2, 'Precipitation should add 1 layer, plus 1 padding layer');
            
            % Verify the top layer is actually fresh snow
            tcase.verifyEqual(round(d_out(1), 2), 150, 'Top layer should be fresh snow density (150)');
        end
        
        function test_melt_event(tcase)
            % Test that high temperature causes melt (via thermo -> melt sub-functions)
            tcase.CF.T_air = 280; % Hot air
            tcase.CF.dsw = 800;   % Intense sun
            tcase.t_vec(1) = 273.15; % Surface at melting point
            
            % Increase timestep to allow significant energy input
            tcase.CF.dt = 3600 * 3; 
            
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, m_tot, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.EC_prev, tcase.M_surf_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyTrue(m_tot > 0, 'High energy input should generate melt');
        end
        
        function test_densification_compaction(tcase)
            % Test that densification returns positive compaction value
            tcase.MP.densification_method = "HerronLangway"; 
            
            % Set density low so there is room to compact
            tcase.d(:) = 300; 
            
            [~, ~, d_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, comp_dens, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.EC_prev, tcase.M_surf_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyTrue(comp_dens > 0, 'Densification should result in positive dry compaction');
            tcase.verifyTrue(all(d_out >= 300), 'Density should increase');
        end
    end
end