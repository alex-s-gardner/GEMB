classdef test_calculate_albedo < matlab.unittest.TestCase
    
    properties
        % Common inputs
        n = 10;
        temperature
        dz
        density
        water
        grain_radius
        albedo_in
        albedo_diffuse_in
        evaporation_condensation
        melt_surface
        
        % Structs
        CF % ClimateForcingStep
        MP % ModelParam
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize standard profiles (Cold Firn Column)
            tcase.temperature              = 260 * ones(tcase.n, 1);  % kelvin
            tcase.dz                       = 0.1 * ones(tcase.n, 1);  % meters
            tcase.density                  = 400 * ones(tcase.n, 1);  % kg m-3
            tcase.water                    = zeros(tcase.n, 1);       % Liquid water content
            tcase.grain_radius             = 0.5 * ones(tcase.n, 1);  % Grain radius [mm]
            tcase.albedo_in                = 0.8 * ones(tcase.n, 1);  % Previous albedo
            tcase.albedo_diffuse_in        = 0.8 * ones(tcase.n, 1);
            tcase.evaporation_condensation = 0;                       % Evap/Cond
            tcase.melt_surface             = 0;                       % Surface Melt
            
            % Initialize ClimateForcingStep Defaults
            tcase.CF.dt                      = 3600; % Time step [s]
            tcase.CF.black_carbon_snow       = 0.1;  % ppmw
            tcase.CF.black_carbon_ice        = 0.1;  % ppmw
            tcase.CF.solar_zenith_angle      = 60;   % degrees
            tcase.CF.cloud_optical_thickness = 1.0; 
            tcase.CF.cloud_fraction          = 0.5;
            tcase.CF.precipitation           = 0;    % Precip
            
            % Initialize ModelParam Defaults
            tcase.MP.albedo_method            = "None";
            tcase.MP.albedo_fixed             = 0.8;
            tcase.MP.albedo_density_threshold = 1023; % Threshold to force fixed albedo
            tcase.MP.density_ice              = 917;
            tcase.MP.albedo_snow              = 0.85;
            tcase.MP.albedo_ice               = 0.48;
            
            % Bougamont2005 specific params
            tcase.MP.albedo_wet_snow_t0 = 15;
            tcase.MP.albedo_dry_snow_t0 = 30;
            tcase.MP.albedo_K           = 7;
            
            % Add source path if needed (Optional)
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
        end
    end
    
    methods (Test)
        
        function test_method_none(tcase)
            % Test: 'None' should return fixed albedo
            tcase.MP.albedo_method = "None";
            tcase.MP.albedo_fixed  = 0.75;
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            tcase.verifyEqual(a_out(1), 0.75, 'AbsTol', 1e-10, 'Method None should return fixed albedo');
        end
        
        function test_density_threshold_override(tcase)
            % Test: High density override should trigger fixed albedo regardless of method
            tcase.MP.albedo_method = "GardnerSharp";
            tcase.MP.albedo_fixed  = 0.4;
            
            % Set threshold low, and density high to trigger override
            tcase.MP.albedo_density_threshold = 300; 
            tcase.density(1) = 350; % Higher than threshold
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            tcase.verifyEqual(a_out(1), 0.4, 'AbsTol', 1e-10, 'Density > Threshold should use fixed albedo');
        end
        
        function test_gardner_sharp(tcase)
            % Test: GardnerSharp functionality (checking range validity)
            tcase.MP.albedo_method = "GardnerSharp";
            
            [a_out, a_diff_out] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            % Basic physical checks
            tcase.verifyTrue(a_out(1) >= 0 && a_out(1) <= 1, 'Albedo must be between 0 and 1');
            tcase.verifyTrue(a_diff_out(1) >= 0 && a_diff_out(1) <= 1, 'Diffuse Albedo must be between 0 and 1');
        end
        
        function test_brune_lefebre(tcase)
            % Test: BruneLeFebre method (spectral bands)
            tcase.MP.albedo_method = "BruneLeFebre";
            tcase.grain_radius(1) = 0.5; % 0.5 mm
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            % We verify it runs and changes the value from input
            tcase.verifyNotEqual(a_out(1), tcase.albedo_in(1), 'BruneLeFebre should update albedo');
            tcase.verifyTrue(a_out(1) > 0.6, 'Fine grain snow should have relatively high albedo');
        end
        
        function test_greuell_konzelmann(tcase)
            % Test: GreuellKonzelmann (Density + Clouds)
            tcase.MP.albedo_method = "GreuellKonzelmann";
            
            % Parameters
            tcase.density(1)        = 450;
            tcase.MP.albedo_ice     = 0.5;
            tcase.MP.albedo_snow    = 0.85;
            tcase.MP.density_ice    = 900;
            tcase.CF.cloud_fraction = 0.5;
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            % Expected Calculation:
            % a = a_ice + (density - density_ice) * (a_snow - a_ice) / (d_fresh - density_ice) + 0.05*(cloud - 0.5)
            % a = 0.5 + (450 - 900) * (0.85 - 0.5) / (300 - 900) + 0
            % a = 0.5 + (-450) * (0.35) / (-600)
            % a = 0.5 + 0.75 * 0.35 = 0.5 + 0.2625 = 0.7625
            expected = 0.5 + (450 - 900) * (0.85 - 0.5) / (300 - 900);
            
            tcase.verifyEqual(a_out(1), expected, 'AbsTol', 1e-5, 'GreuellKonzelmann calc mismatch');
        end
        
        function test_bougamont_bamber_decay(tcase)
            % Test: Bougamont2005 time decay
            tcase.MP.albedo_method = "Bougamont2005";
            
            % Cold snow (temperature < -10C), Dry (water=0)
            tcase.temperature(1) = 250; % -23 C
            tcase.water(1) = 0;
            tcase.albedo_in(1) = 0.8;
            tcase.CF.dt = 86400; % 1 day
            
            % Constants
            % albedo_K = 7
            % dry_snow_t0 = 30
            % t0 = 10 * K + dry_snow_t0 = 70 + 30 = 100 days
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            t0 = 100;
            da = (tcase.albedo_in(1) - tcase.MP.albedo_ice) / t0 * 1; % dt=1 day
            expected = tcase.albedo_in(1) - da;
            
            tcase.verifyEqual(a_out(1), expected, 'AbsTol', 1e-5, 'Bougamont2005 decay incorrect');
        end
        
        function test_bougamont_bamber_snowfall(tcase)
            % Test: Bougamont2005 albedo reset due to snowfall (via Condensation logic in albedo.m)
            tcase.MP.albedo_method         = "Bougamont2005";
            tcase.albedo_in(1)             = 0.6; % Old snow
            tcase.evaporation_condensation = 5;   % Significant condensation/deposition
            tcase.temperature(1)           = 260; % Cold enough for solid deposition
            tcase.CF.precipitation         = 0;   % No precip initially
            
            % Logic: precipitation becomes precipitation + (evaporation_condensation/300)*1000
            % a = a_snow - (a_snow - a_old) * exp(-precipitation/z_snow)
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            tcase.verifyTrue(a_out(1) > tcase.albedo_in(1), 'Deposition should increase albedo');
        end
        
        function test_thin_layer_mixing(tcase)
            % Test: GardnerSharp with thin snow layer over ice (Alexander 2014 logic)
            tcase.MP.albedo_method = "GardnerSharp";
            
            % Setup thin snow layer
            tcase.dz(1)      = 0.05;  % 5cm ( < 10cm threshold)
            tcase.density(1) = 350;   % Snow
            tcase.density(2) = 900;   % Ice (assume pore closeoff density is 830)
            
            [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
            
            % Just verify that mixing logic ran (result should not be pure snow albedo)
            % Validating exact number is hard without reproducing `albedo_gardner`,
            % but we can check it returns a valid number.
            tcase.verifyTrue(a_out(1) >= 0 && a_out(1) <= 1);
        end
        
        function test_ice_albedo_aging(tcase)
             % Test: Ice surface logic (density > density_ice)
             % Logic: a = max(min_ice + (max_ice - min_ice)*exp(-M/200), min_ice)
             tcase.MP.albedo_method = "GardnerSharp";
             tcase.density(1) = 917; % Ice
             
             tcase.melt_surface = 100; % 100 kg m-2 melt
             tcase.water(1) = 0;
             
             % Constants from albedo.m
             % ice_max = 0.58
             % ice_min = MP.albedo_ice (0.48)
             % K = 200
             
             [a_out, ~] = calculate_albedo(tcase.temperature, tcase.dz, tcase.density, tcase.water, tcase.grain_radius, ...
                tcase.albedo_in, tcase.albedo_diffuse_in, tcase.evaporation_condensation, tcase.melt_surface, tcase.CF, tcase.MP);
             
             M = 100;
             expected = 0.48 + (0.58 - 0.48) * exp(-1.0 * (M/200.0));
             
             tcase.verifyEqual(a_out(1), expected, 'AbsTol', 1e-5, 'Ice albedo aging incorrect');
        end
    end
end