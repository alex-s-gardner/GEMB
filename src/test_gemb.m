classdef test_gemb < matlab.unittest.TestCase
    
    properties
        % Input Structures
        cf % ClimateForcing
        mp % ModelParam
        
        % Column State
        t_vec
        dz
        d
        w
        re
        gdn
        gsp
        a
        a_diff
        
        % Time settings
        dt_days = 1/24; % 1 hour
        n_steps = 5;
    end
    
    methods (TestClassSetup)
        function setup_environment(tcase)
            % 1. Create Mock for 'model_initialize_output.m'
            fid = fopen('model_initialize_output.m', 'w');
            fprintf(fid, 'function [idx, out_data, out_cum] = model_initialize_output(len, cf, mp)\n');
            fprintf(fid, '    %% Mock output initializer\n');
            fprintf(fid, '    idx = true(length(cf.daten), 1);\n');
            fprintf(fid, '    out_cum.count = 0;\n');
            fprintf(fid, '    %% Variables defined in gemb.m for eval()\n');
            fprintf(fid, '    vars = {"sw_net","lw_net","shf","lhf","M","R","F","EC","Ra","M_added","compaction_dens","compaction_melt","d1","a1","re1","Q_net","FAC"};\n');
            fprintf(fid, '    for i=1:length(vars), out_cum.(vars{i}) = 0; end\n');
            fprintf(fid, '    \n');
            fprintf(fid, '    %% Initialize Output Struct\n');
            fprintf(fid, '    n_time = length(cf.daten);\n');
            fprintf(fid, '    padding = mp.output_padding + len;\n');
            fprintf(fid, '    out_data.time = cf.daten;\n');
            fprintf(fid, '    out_data.dz = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.d = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.T = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.W = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.re = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.gdn = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.gsp = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.ps = zeros(padding, n_time);\n');
            fprintf(fid, '    out_data.m = zeros(n_time, 1);\n');
            fprintf(fid, '    for i=1:length(vars), out_data.(vars{i}) = zeros(n_time, 1); end\n');
            fprintf(fid, 'end\n');
            fclose(fid);
            
            % 2. Create Mock for 'gemb_core.m'
            fid = fopen('gemb_core.m', 'w');
            fprintf(fid, 'function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, sw_net, shf, lhf, ulw, Ra, M, R, F, M_added, E_added, compaction_dens, compaction_melt] = gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, ClimateForcingStep, ModelParam, verbose)\n');
            fprintf(fid, '    %% Mock Physics Kernel\n');
            fprintf(fid, '    sw_net = 0; shf = 0; lhf = 0; \n');
            fprintf(fid, '    ulw = ClimateForcingStep.dlw; %% Balance LW to prevent T drift error\n');
            fprintf(fid, '    Ra = 0; M = 0; R = 0; F = 0; M_added = 0; E_added = 0;\n');
            fprintf(fid, '    compaction_dens = 0; compaction_melt = 0;\n');
            fprintf(fid, 'end\n');
            fclose(fid);

            % 3. Setup Path
            addpath(pwd);
            import matlab.unittest.fixtures.PathFixture
            try
                tcase.applyFixture(PathFixture('../src'));
            catch
                addpath('src');
            end
            
            tcase.addTeardown(@() delete('model_initialize_output.m'));
            tcase.addTeardown(@() delete('gemb_core.m'));
        end
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % 1. Create Mock ClimateForcing (Time Series)
            tcase.cf.daten = (0:tcase.n_steps-1)' * tcase.dt_days;
            
            n = tcase.n_steps;
            tcase.cf.T_air0 = 260 * ones(n, 1);
            tcase.cf.V0 = 5 * ones(n, 1);
            tcase.cf.dlw0 = 200 * ones(n, 1);
            tcase.cf.dsw0 = 0 * ones(n, 1); 
            tcase.cf.e_air0 = 100 * ones(n, 1);
            tcase.cf.p_air0 = 100000 * ones(n, 1);
            tcase.cf.P0 = zeros(n, 1); 
            
            % Location Params
            tcase.cf.Vz = 2;
            tcase.cf.Tz = 2;
            tcase.cf.T_air_mean = 260;
            tcase.cf.P_mean = 0;
            tcase.cf.V_mean = 5;
            tcase.cf.elev = 2000;
            tcase.cf.lat = 70;
            tcase.cf.lon = -40;
            
            % 2. Create Mock ModelParam
            tcase.mp.run_prefix = "TEST";
            tcase.mp.n_spinup_cycles = 0;
            tcase.mp.density_ice = 917;
            tcase.mp.output_padding = 10;
            
            % Radiative Params
            tcase.mp.black_carbon_snow = 0;
            tcase.mp.black_carbon_ice = 0;
            tcase.mp.cloud_optical_thickness = 0;
            tcase.mp.solar_zenith_angle = 0;
            tcase.mp.dsw_diffuse = 0;
            tcase.mp.cloud_fraction = 0;
            
            % 3. Initialize Column State
            cols = 10;
            tcase.t_vec = 260 * ones(cols, 1);
            tcase.dz = 0.1 * ones(cols, 1);
            tcase.d = 400 * ones(cols, 1);
            tcase.w = zeros(cols, 1);
            tcase.re = 0.5 * ones(cols, 1);
            tcase.gdn = 0.5 * ones(cols, 1);
            tcase.gsp = 0.5 * ones(cols, 1);
            tcase.a = 0.8 * ones(cols, 1);
            tcase.a_diff = 0.8 * ones(cols, 1);
        end
    end
    
    methods (Test)
        
        function test_basic_execution(tcase)
            % Verify the model runs through 5 timesteps without crashing
            out_data = gemb(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.cf, tcase.mp, false);
            
            tcase.verifyTrue(isstruct(out_data));
            
            % Fix: Check the number of columns (time steps), not length (max dim)
            % OutData.T is [layers x time]
            tcase.verifyEqual(size(out_data.T, 2), tcase.n_steps, ...
                'Output time dimension should match number of steps');
        end
        
        function test_spinup_logic(tcase)
            % Test 1 spinup cycle
            tcase.mp.n_spinup_cycles = 1;
            
            evalc('out_data = gemb(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, tcase.cf, tcase.mp, false);');
            
            % Ensure output structure contains only the final cycle
            % out_data.m is a vector [1 x time]
            tcase.verifyEqual(length(out_data.m), tcase.n_steps, ...
                'Output should only contain data for the final simulation cycle');
        end
        
        function test_mass_conservation_check(tcase)
            % gemb.m throws error if M_change ~= 0
            try
                gemb(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                    tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                    tcase.cf, tcase.mp, false);
            catch ME
                tcase.verifyFail(['Mass conservation check failed: ' ME.message]);
            end
        end
        
        function test_bottom_temperature_check(tcase)
            % Verify error if bottom temperature changes
            fid = fopen('gemb_core.m', 'w');
            fprintf(fid, 'function [T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, sw_net, shf, lhf, ulw, Ra, M, R, F, M_added, E_added, compaction_dens, compaction_melt] = gemb_core(T, dz, d, W, re, gdn, gsp, a, a_diffuse, EC, M_surf, ClimateForcingStep, ModelParam, verbose)\n');
            fprintf(fid, '    T(end) = T(end) + 10; %% Force Change\n'); 
            fprintf(fid, '    sw_net=0; shf=0; lhf=0; ulw=ClimateForcingStep.dlw; Ra=0; M=0; R=0; F=0; M_added=0; E_added=0; compaction_dens=0; compaction_melt=0;\n');
            fprintf(fid, 'end\n');
            fclose(fid);
            
            verifyError(tcase, @() gemb(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diff, ...
                tcase.cf, tcase.mp, false), ...
                ?MException);
        end
    end
end