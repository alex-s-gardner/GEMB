classdef test_calculate_accumulation < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 5;
        t_vec
        dz
        density
        water
        grain_radius
        grain_dendricity
        grain_sphericity
        albedo_in
        albedo_diffuse_in
        
        % Constants and Parameters
        density_ice = 917;
        temperature_air_mean = 260;
        precipitation_mean = 200;
        wind_speed_mean = 5;
        alb_snow = 0.85;
        alb_method = "GardnerSharp"; % Updated to string
        dz_min = 0.05; % Minimum thickness to create new layer
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)

            % Initialize standard profiles (cold firn column)
            tcase.t_vec             = 260 * ones(tcase.n, 1);
            tcase.dz                = 0.1 * ones(tcase.n, 1);
            tcase.density           = 400 * ones(tcase.n, 1); 
            tcase.water             = zeros(tcase.n, 1);
            tcase.grain_radius      = 0.5 * ones(tcase.n, 1);
            tcase.grain_dendricity  = 0.5 * ones(tcase.n, 1);
            tcase.grain_sphericity  = 0.5 * ones(tcase.n, 1);
            tcase.albedo_in         = 0.7 * ones(tcase.n, 1);
            tcase.albedo_diffuse_in = 0.7 * ones(tcase.n, 1);
            
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
            
            % Setup Structures
            CF.precipitation        = 0;
            CF.temperature_air      = 270;
            CF.wind_speed           = 5;
            CF.precipitation_mean   = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean      = tcase.wind_speed_mean;
            
            MP.new_snow_method            = "150kgm2";
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15 ;
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, ra_out] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            
            tcase.verifyEqual(dz_out, tcase.dz);
            tcase.verifyEqual(d_out, tcase.density);
            tcase.verifyEqual(t_out, tcase.t_vec);
            tcase.verifyEqual(ra_out, 0);
        end
        
        function test_large_snow_event_new_layer(tcase)
            % Large snowfall (> dz_min) should create a new layer on top
            density_snow = 150;
            
            % Setup Structures
            CF.precipitation        = 50; % kg m-2
            CF.temperature_air      = 260; % Cold air -> Snow
            CF.wind_speed           = 5;
            CF.precipitation_mean   = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean      = tcase.wind_speed_mean;
            
            MP.new_snow_method            = "150kgm2";
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15 ;
            
            [t_out, dz_out, d_out, ~, ~, gdn_out, gsp_out, a_out, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            
            expected_dz = CF.precipitation / density_snow;

            % Check vector growth
            tcase.verifyEqual(length(dz_out), tcase.n + 1, 'Large snow should add a layer');
            
            % Check top layer properties
            tcase.verifyEqual(d_out(1), density_snow, 'Top layer should have fresh snow density');

            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10);
            tcase.verifyEqual(t_out(1), CF.temperature_air);
            tcase.verifyEqual(a_out(1), tcase.alb_snow);
            
            % Check default microstructure for new snow
            tcase.verifyEqual(gdn_out(1), 1.0); % Dendricity default
            tcase.verifyEqual(gsp_out(1), 0.5); % Sphericity default
        end
        
        function test_small_snow_event_merge(tcase)

            % Small snowfall (< dz_min) should merge with top layer
            density_snow = 150;
            
            % Setup Structures
            CF.precipitation        = 2; % 2 kg m-2
            CF.temperature_air      = 260;
            CF.wind_speed           = 5;
            CF.precipitation_mean   = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean      = tcase.wind_speed_mean;
            
            MP.new_snow_method            = "150kgm2";
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15 ;
            
            old_mass = tcase.density(1) * tcase.dz(1);
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, a_out, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            
            % Verify no new layer
            tcase.verifyEqual(length(dz_out), tcase.n, 'Small snow should merge');
            
            % Verify mass conservation and mixing
            new_mass    = old_mass + CF.precipitation;
            expected_dz = tcase.dz(1) + CF.precipitation/density_snow;
            expected_d  = new_mass / expected_dz;
            
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Depth should increase');
            tcase.verifyEqual(d_out(1),   expected_d, 'AbsTol', 1e-10, 'Density should decrease (mix with light snow)');
            
            % Verify Temperature weighting
            expected_t = (CF.temperature_air * CF.precipitation + tcase.t_vec(1) * old_mass) / new_mass;
            tcase.verifyEqual(t_out(1), expected_t, 'AbsTol', 1e-10, 'Temperature should be weighted average');
            
            % Verify Albedo update
            expected_a = (tcase.alb_snow * CF.precipitation + tcase.albedo_in(1) * old_mass) / new_mass;
            tcase.verifyEqual(a_out(1), expected_a, 'AbsTol', 1e-10, 'Albedo should be weighted average');
        end
        
        function test_rain_event(tcase)
            % Rain (temperature_air > 0C) should add mass and heat to top layer, not thickness
            
            % Setup Structures
            CF.precipitation        = 10;
            CF.temperature_air      = 275; % > 273.15 -> Rain
            CF.wind_speed           = 0;
            CF.precipitation_mean   = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean      = tcase.wind_speed_mean;
            
            MP.new_snow_method            = "150kgm2";
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15;
            
            % Constants from function
            lf = 0.3345e6;
            ci = 2102;
            
            old_mass = tcase.density(1) * tcase.dz(1);
            
            [t_out, dz_out, d_out, ~, ~, ~, ~, ~, ~, ra_out] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            
            % Rain Output flag
            tcase.verifyEqual(ra_out, CF.precipitation);
            
            % Mass update (thickness stays same, density increases for rain in this model)
            % The code: dz(1) stays same, density(1) = mass/dz(1)
            new_mass = old_mass + CF.precipitation;
            tcase.verifyEqual(d_out(1), new_mass/tcase.dz(1), 'AbsTol', 1e-10);
            tcase.verifyEqual(dz_out(1), tcase.dz(1), 'Thickness unchanged for rain unless ice density hit');
            
            % Temperature update (includes Latent Heat logic)
            % T(1) = (precipitation *(temperature_air + LF/CI) + T(1) * mInit(1)) / mass;
            term_rain  = CF.precipitation * (CF.temperature_air + lf/ci);
            term_snow  = tcase.t_vec(1) * old_mass;
            expected_t = (term_rain + term_snow) / new_mass;
            
            tcase.verifyEqual(t_out(1), expected_t, 'AbsTol', 1e-8, 'Temperature should account for rain latent heat');
        end
        
        function test_rain_density_cap(tcase)
            % If rain causes density to exceed ice density, it should clamp density and increase thickness
            
            % Setup Structures
            CF.precipitation        = 500; % Huge rain event
            CF.temperature_air      = 275;
            CF.wind_speed           = 0;
            CF.precipitation_mean   = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean      = tcase.wind_speed_mean;
            
            MP.new_snow_method            = "150kgm2";
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15 ;
            
            % Start with density near ice
            tcase.density(1) = 900; 
            tcase.dz(1) = 0.1; % Mass = 90
            
            % New mass = 590. If dz kept 0.1, density = 5900 (impossible).
            % Expect density -> density_ice, dz -> increases.
            
            [~, dz_out, d_out, ~, ~, ~, ~, ~, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            
            tcase.verifyEqual(d_out(1), tcase.density_ice, 'AbsTol', 1e-10, 'Density should be capped at ice density');
            
            total_mass  = 90 + 500;
            expected_dz = total_mass / tcase.density_ice;
            tcase.verifyEqual(dz_out(1), expected_dz, 'AbsTol', 1e-10, 'Thickness should expand to conserve mass at ice density');
        end
        
        function test_density_methods(tcase)
            % Test the different switches for new snow density
            
            % Common setup
            CF.precipitation = 50; % Ensure new layer
            CF.temperature_air = 250;
            CF.wind_speed = 5;
            CF.precipitation_mean = tcase.precipitation_mean;
            CF.temperature_air_mean = tcase.temperature_air_mean;
            CF.wind_speed_mean = tcase.wind_speed_mean;
            
            MP.column_dzmin               = tcase.dz_min;
            MP.density_ice                = tcase.density_ice;
            MP.albedo_snow                = tcase.alb_snow;
            MP.albedo_method              = tcase.alb_method;
            MP.rain_temperature_threshold = 273.15 ;
            
            % Method 1: Antarctica (350) -> "350kgm2"
            MP.new_snow_method = "350kgm2";
            [~, ~, d1, ~, ~, ~, ~, ~, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            tcase.verifyEqual(d1(1), 350, 'Method "350kgm2" should be 350');
            
            % Method 2: Greenland -> "Fausto"
            MP.new_snow_method = "Fausto";
            [~, ~, d2, ~, ~, ~, ~, ~, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            tcase.verifyEqual(d2(1), 315, 'Method "Fausto" should be 315');
            
            % Method 3: Kaspers -> "Kaspers"
            % dSnow=(7.36e-2 + 1.06e-3*min(temperature_air_mean,273.15) + 6.69e-2*precipitation_mean/1000. + 4.77e-3*wind_speed_mean)*1000.
            MP.new_snow_method = "Kaspers";
            expected_3 = (7.36e-2 + 1.06e-3*tcase.temperature_air_mean + 6.69e-2*tcase.precipitation_mean/1000 + 4.77e-3*tcase.wind_speed_mean)*1000;
            [~, ~, d3, ~, ~, ~, ~, ~, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            tcase.verifyEqual(d3(1), expected_3, 'AbsTol', 1e-5, 'Method "Kaspers" calculation incorrect');
            
            % Method 4: Kuipers Munneke -> "KuipersMunneke"
            % dSnow = 481.0 + 4.834*(temperature_air_mean-273.15);
            MP.new_snow_method = "KuipersMunneke";
            expected_4 = 481.0 + 4.834*(tcase.temperature_air_mean - 273.15);
            [~, ~, d4, ~, ~, ~, ~, ~, ~, ~] = calculate_accumulation(...
                tcase.t_vec, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.grain_dendricity, tcase.grain_sphericity, tcase.albedo_in, tcase.albedo_diffuse_in, ...
                CF, MP, false);
            tcase.verifyEqual(d4(1), expected_4, 'AbsTol', 1e-5, 'Method "KuipersMunneke" calculation incorrect');
        end
    end
end