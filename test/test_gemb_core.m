classdef test_gemb_core < matlab.unittest.TestCase
    
    properties
        % Grid state vectors
        n = 10;
        t_vec
        dz
        density
        water
        grain_radius
        grain_density
        grain_sphericity
        albedo
        albedo_diffuse
        
        % State scalars
        evaporation_condensation_prev = 0;
        mass_surface_prev = 0;
        
        % Structs
        CF % ClimateForcingStep
        MP % ModelParam
        
        % CHANGED: Set to true to avoid scope error in gemb_core.m 
        % (CI is defined only when verbose=true)
        verbose = true; 
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
            tcase.dz               = 0.08 * ones(tcase.n, 1);
            tcase.density          = 400 * ones(tcase.n, 1);
            tcase.water            = zeros(tcase.n, 1);
            tcase.grain_radius     = 0.5 * ones(tcase.n, 1);
            tcase.grain_density    = 0.5 * ones(tcase.n, 1);
            tcase.grain_sphericity = 0.5 * ones(tcase.n, 1);
            tcase.albedo           = 0.8 * ones(tcase.n, 1);
            tcase.albedo_diffuse   = 0.8 * ones(tcase.n, 1);
            
            % --- Initialize ClimateForcingStep (CF) ---
            tcase.CF.dt                      = 3600; % 1 hour
            tcase.CF.precipitation           = 0;
            tcase.CF.temperature_air         = 265;
            tcase.CF.wind_speed              = 5;
            tcase.CF.vapor_pressure          = 400;
            tcase.CF.pressure_air            = 100000;
            tcase.CF.longwave_downward       = 300;
            tcase.CF.shortwave_downward      = 200;
            tcase.CF.shortwave_diffuse_diff  = 50;
            tcase.CF.black_carbon_snow       = 0;
            tcase.CF.black_carbon_ice        = 0;
            tcase.CF.solar_zenith_angle      = 60;
            tcase.CF.cloud_optical_thickness = 0;
            tcase.CF.cloud_fraction          = 0;
            
            % Location/Mean parameters
            tcase.CF.wind_observation_height        = 2;
            tcase.CF.temperature_observation_height = 2;
            tcase.CF.temperature_air_mean           = 260;
            tcase.CF.precipitation_mean             = 200;
            tcase.CF.wind_speed_mean                = 5;
            
            % --- Initialize ModelParam (MP) ---
            tcase.MP.albedo_method = "GardnerSharp"; 
            tcase.MP.albedo_ice    = 0.45;
            tcase.MP.albedo_snow   = 0.85;
            tcase.MP.albedo_fixed  = 0.7;
            
            % "desnity" typo matches source code expectation in albedo.m
            tcase.MP.albedo_density_threshold = 1023; 
            
            tcase.MP.albedo_wet_snow_t0          = 15;
            tcase.MP.albedo_dry_snow_t0          = 30;
            tcase.MP.albedo_K                    = 7;
            tcase.MP.shortwave_absorption_method = 1;
            
            % Emissivity
            tcase.MP.emissivity                      = 0.98;
            tcase.MP.emissivity_grain_radius_large     = 0.97; % Required by thermo.m
            tcase.MP.emissivity_method               = "uniform"; % String required by thermo.m
            tcase.MP.longwave_upward_delta           = 0;
            tcase.MP.emissivity_grain_radius_threshold = 10;
            
            % Thermal & Grid
            tcase.MP.thermal_conductivity_method  = "Sturm"; 
            tcase.MP.column_dzmin                 = 0.05;
            tcase.MP.column_dzmax                 = 0.10; 
            tcase.MP.new_snow_method              = "150kgm2"; 
            tcase.MP.density_ice                  = 917;
            tcase.MP.water_irreducible_saturation = 0.07;
            
            % Densification
            tcase.MP.densification_method     = "HerronLangway"; 
            tcase.MP.densification_coeffs_M01 = "Gre_RACMO_GS_SW0"; % Default needed
            
            % Set zmax slightly larger than initial depth to prevent removal.
            % Initial depth = 0.8m. 
            tcase.MP.column_zmax = 0.9; 
            tcase.MP.column_zmin = 0.5; 
            tcase.MP.column_ztop = 2;
            tcase.MP.column_zy   = 1.1;

            tcase.MP.surface_roughness_effective_ratio = 0.1; 
            tcase.MP.rain_temperature_threshold        = 273.15; 


            tcase.MP.surface_roughness_effective_ratio = 0.1; 
            tcase.MP.rain_temperature_threshold = 273.15; 

            % --- CRITICAL FIX: Add dt_divisors ---
            % thermo.m requires this field to determine stable time steps.
            % Normally generated by gemb.m using fast_divisors.
            % Hardcoding common divisors for 3600s (1 hour)
            tcase.MP.dt_divisors = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, ...
                                    24, 25, 30, 32, 36, 40, 45, 48, 50, 60, 72, 75, ...
                                    80, 90, 100, 120, 144, 150, 180, 200, 225, 240, ...
                                    300, 360, 400, 450, 600, 720, 900, 1200, 1800, 3600];
        end
    end
    
    methods (Test)
        
        function test_pipeline_execution(tcase)
            % Basic "Smoke Test": Verify it runs without error
            % Maps to 23 outputs in gemb_core
            
            [t_out, dz_out, density_out, ~, ~, ~, ~, a_out, ~, ~, ~, shortwave_net, ...
             heat_flux_sensible, heat_flux_latent, longwave_upward, ~, m_tot, r_tot, f_tot, m_add, e_add, comp_dens, comp_melt] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_density, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, tcase.evaporation_condensation_prev, tcase.mass_surface_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            % Check Outputs
            % Size logic: 
            % Initial 10 layers (0.8m). zmax=0.9m.
            % layer_management sees depth < zmax -> adds 1 layer.
            % Result: 11 layers.
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Expect 1 padding layer');
            
            tcase.verifyEqual(size(dz_out), size(t_out));
            tcase.verifyEqual(size(dz_out), size(density_out));
            
            % Fluxes
            tcase.verifyTrue(isscalar(shortwave_net));
            tcase.verifyTrue(isscalar(heat_flux_sensible));
            tcase.verifyTrue(isscalar(heat_flux_latent));
            tcase.verifyTrue(isscalar(longwave_upward));
            
            % Verify compaction variables are returned
            tcase.verifyTrue(isscalar(comp_dens));
            tcase.verifyTrue(isscalar(comp_melt));
        end
        
        function test_accumulation_event(tcase)
            % Test that precipitation adds a layer
            
            % 1. Use precipitation = 10 kg/m2.
            tcase.CF.precipitation   = 10; 
            tcase.CF.temperature_air = 260; % Snow
            
            % 2. Set zmax huge to force padding logic instead of removal logic.
            tcase.MP.column_zmax = 100.0;
            
            % Capture 23 outputs (using ~ for unused ones)
            [~, dz_out, density_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_density, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, tcase.evaporation_condensation_prev, tcase.mass_surface_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            % Expected Layers:
            % 10 (Original) + 1 (Snow Accumulation) + 1 (Zmax Padding) = 12
            tcase.verifyEqual(length(dz_out), tcase.n + 2, 'Precipitation should add 1 layer, plus 1 padding layer');
            
            % Verify the top layer is actually fresh snow
            tcase.verifyEqual(round(density_out(1), 2), 150, 'Top layer should be fresh snow density (150)');
        end
        
        function test_melt_event(tcase)
            % Test that high temperature causes melt (via thermo -> melt sub-functions)
            tcase.CF.temperature_air    = 280; % Hot air
            tcase.CF.shortwave_downward = 800;   % Intense sun
            tcase.t_vec(1) = 273.15; % Surface at melting point
            
            % Increase timestep to allow significant energy input (10800s = 3 hours)
            % Note: 3600 is a divisor of 10800, so our hardcoded dt_divisors list works.
            tcase.CF.dt = 3600 * 3; 
            
            % Capture m_tot (index 17)
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, m_tot, ~, ~, ~, ~, ~, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_density, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, tcase.evaporation_condensation_prev, tcase.mass_surface_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyTrue(m_tot > 0, 'High energy input should generate melt');
        end
        
        function test_densification_compaction(tcase)
            % Test that densification returns positive compaction value
            tcase.MP.densification_method = "HerronLangway"; 
            
            % Set density low so there is room to compact
            tcase.density(:) = 300; 
            
            % Capture density_out (3) and comp_dens (22)
            [~, ~, density_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, comp_dens, ~] = ...
                gemb_core(tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_density, tcase.grain_sphericity, tcase.albedo, tcase.albedo_diffuse, tcase.evaporation_condensation_prev, tcase.mass_surface_prev, ...
                tcase.CF, tcase.MP, tcase.verbose);
            
            tcase.verifyTrue(comp_dens > 0, 'Densification should result in positive dry compaction');
            tcase.verifyTrue(all(density_out >= 300), 'Density should increase');
        end
    end
end