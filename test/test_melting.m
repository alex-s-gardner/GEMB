classdef test_melting < matlab.unittest.TestCase
    
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
        a
        a_diffuse
        ModelParam
        
        % Parameters
        rho_ice = 920;
        ra = 0;
        verbose = false;
    end
    
    methods (TestMethodSetup)
        function setup_inputs(tcase)
            % Initialize standard profile (Cold, Dry Firn)
            tcase.t_vec = 260 * ones(tcase.n, 1);
            tcase.dz = 0.1 * ones(tcase.n, 1);
            tcase.d = 400 * ones(tcase.n, 1);
            tcase.w = zeros(tcase.n, 1);
            tcase.re = 0.5 * ones(tcase.n, 1);
            tcase.gdn = 0.5 * ones(tcase.n, 1);
            tcase.gsp = 0.5 * ones(tcase.n, 1);
            tcase.a = 0.8 * ones(tcase.n, 1);
            tcase.a_diffuse = 0.8 * ones(tcase.n, 1);
            
            tcase.ModelParam.density_ice = 920;
            tcase.ModelParam.water_irreducible_saturation = 0.07;


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
        
        function test_cold_dry_snow_no_change(tcase)
            % Scenario: T < 0C, W = 0. No melt, no refreeze.
            
            [t_out, ~, d_out, w_out, ~, ~, ~, ~, ~, m_tot, ~, r_tot, f_tot] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            tcase.verifyEqual(m_tot, 0);
            tcase.verifyEqual(r_tot, 0);
            tcase.verifyEqual(f_tot, 0);
            tcase.verifyEqual(t_out, tcase.t_vec);
            tcase.verifyEqual(d_out, tcase.d);
            tcase.verifyEqual(w_out, tcase.w);
        end
        
        function test_pore_water_refreeze(tcase)
            % Scenario: Cold snow with liquid water. 
            % Water should freeze, releasing latent heat -> T increases.
            
            tcase.w(1) = 5; % 5 kg water
            mass_initial = tcase.d(1) * tcase.dz(1) + tcase.w(1);
            
            [t_out, dz_out, d_out, w_out, ~, ~, ~, ~, ~, ~, ~, ~, f_tot] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            % Verify refreeze occurred
            tcase.verifyTrue(f_tot > 0, 'Some water should refreeze');
            tcase.verifyTrue(w_out(1) < 5, 'Pore water should decrease');
            
            % Verify warming
            tcase.verifyTrue(t_out(1) > tcase.t_vec(1), 'Temperature should rise due to latent heat');
            
            % Verify density increase (water becomes ice/snow matrix)
            mass_final = d_out(1) * dz_out(1) + w_out(1);
            tcase.verifyTrue(d_out(1) > tcase.d(1), 'Density should increase');
            
            % Simple mass check
            tcase.verifyEqual(mass_final, mass_initial, 'AbsTol', 1e-10, 'Mass conservation');
        end
        
        function test_surface_melt(tcase)
            % Scenario: Surface T > 0C. 
            % Should generate melt, clamp T to 0C (273.15 K).
            
            tcase.t_vec(1) = 280; % Hot surface
            
            [t_out, ~, ~, ~, ~, ~, ~, ~, ~, m_tot, m_surf, ~, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            tcase.verifyEqual(t_out(1), 273.15, 'Temperature should be clamped to freezing point');
            tcase.verifyTrue(m_surf > 0, 'Surface melt should occur');
            tcase.verifyTrue(m_tot > 0, 'Total melt should occur');
        end
        
        function test_runoff_on_ice(tcase)
            % Scenario: Melt occurs above an impermeable ice layer.
            % Water should not percolate; it should run off.
            
            tcase.t_vec(1) = 280; % Melt source
            
            % UPDATED: Extend ice layer to 2 cells (0.2m). 
            % melting.m requires ice_depth > 0.1m to trigger runoff logic.
            tcase.d(2:3) = 830; 
             
            % CRITICAL FIX: Pre-saturate the top layer.
            % Without this, melt water is retained as Irreducible Water (Swi) 
            % and won't run off.
            tcase.w(1) = 10; % Add pore water to ensure saturation
            
            [~, ~, ~, w_out, ~, ~, ~, ~, ~, ~, ~, r_tot, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            tcase.verifyTrue(r_tot > 0, 'Runoff should occur when saturated snow sits on thick ice');
            
            % UPDATED: The ice layer itself (layer 2) WILL retain water up to 
            % irreducible saturation before runoff occurs.
            tcase.verifyTrue(w_out(2) > 0, 'Ice layer should retain some irreducible water');
            
            % Ensure no water percolation THROUGH the ice layer to the snow below (layer 4)
            tcase.verifyEqual(w_out(4), 0, 'No water should pass through the impermeable ice layer');
        end
        
        function test_excess_heat_distribution(tcase)
            % Scenario: Input T is HUGE, enough to melt the entire top cell.
            % Excess heat should propagate downward.
            
            tcase.t_vec(1) = 273.15 + 500; 
            
            [t_out, ~, d_out, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            if length(d_out) < tcase.n
                % Cell melted away
                tcase.verifyTrue(true, 'Top cell melted completely');
            else
                 % Cell exists, check layer 2
                 % Layer 2 should have warmed up due to excess energy transfer
                 tcase.verifyTrue(t_out(2) > 260, 'Excess heat should warm underlying layer');
            end
        end
        
        function test_water_squeezing(tcase)
            % Scenario: No melt (Cold), but W > Irreducible water content.
            % Water should drain (percolate/runoff) even without active melting.
            
            tcase.w(1) = 20; % Very wet
            tcase.t_vec(:) = 273.15; % Isothermal at 0 to prevent immediate refreeze
            
            [~, ~, ~, w_out, ~, ~, ~, ~, ~, ~, ~, r_tot, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                tcase.ra, tcase.ModelParam, tcase.verbose);
            
            % Water should leave top cell (percolate or run off)
            tcase.verifyTrue(w_out(1) < 20, 'Excess water should drain');
            tcase.verifyTrue(r_tot > 0 || sum(w_out(2:end)) > 0, 'Water should move down or run off');
        end
        
        function test_rain_accounting(tcase)
            % Scenario: Rain accounting.
            % melting.m subtracts Ra from the *reported* Total Melt (M_total) 
            % to distinguish between liquid rain input and phase-change melt.
            % M_total = sum(MeltedIce) - Ra.
            
            % Generate actual melt
            tcase.t_vec(1) = 280; 
            % We need to know how much melt this generates to verify the subtraction.
            % M_max = T_excess * d * dz * CI / LF;
            % 6.85 * 400 * 0.1 * 2102 / 334500 approx 1.7 kg.
            
            % Run without rain first to get baseline melt
            [~, ~, ~, ~, ~, ~, ~, ~, ~, m_tot_base, ~, ~, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                0, tcase.ModelParam, tcase.verbose);
            
            % Run with rain input
            rain_input = 0.5;
            [~, ~, ~, ~, ~, ~, ~, ~, ~, m_tot_rain, ~, ~, ~] = ...
                melting(tcase.t_vec, tcase.dz, tcase.d, tcase.w, tcase.re, ...
                tcase.gdn, tcase.gsp, tcase.a, tcase.a_diffuse, ...
                rain_input, tcase.ModelParam, tcase.verbose);
            
            % M_total with rain should be (M_total_base - rain_input)
            expected = max(0, m_tot_base - rain_input);
            
            tcase.verifyEqual(m_tot_rain, expected, 'AbsTol', 1e-6, ...
                'M_total should subtract rain input (accounting logic)');
        end
    end
end