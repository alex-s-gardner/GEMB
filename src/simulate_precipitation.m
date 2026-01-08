function precip = simulate_precipitation(dec_year, coeffs)
% SIMULATE_PRECIPITATION Generates synthetic rainfall using fitted coeffs.
%
%   INPUTS:
%       dec_year - N x 1 vector of desired decimal years to simulate
%       coeffs   - Struct from fit_precipitation
%
%   OUTPUT:
%       precip - N x 1 vector of simulated precipitation

    dec_year = dec_year(:);
    
    n = length(dec_year);
    precip = zeros(n, 1);
    
    % 1. Reconstruct Time-Varying Parameters
    t_season = mod(dec_year, 1);
    
    % Basis functions
    X = [ones(n,1), sin(2*pi*t_season), cos(2*pi*t_season)];
    
    % Calculate params and clamp to valid ranges
    P01_t = X * coeffs.P01_harmonics;
    P11_t = X * coeffs.P11_harmonics;
    Alpha_t = X * coeffs.Alpha_harmonics;
    Beta_t  = X * coeffs.Beta_harmonics;
    
    % Clamping ensures stability if harmonics swing too wide
    P01_t = max(0, min(1, P01_t));
    P11_t = max(0, min(1, P11_t));
    Alpha_t = max(0.1, Alpha_t); % Shape must be > 0
    Beta_t  = max(0.1, Beta_t);  % Scale must be > 0
    
    % 2. Simulation Loop
    % We must loop sequentially because state depends on previous step
    
    % Initial state (random guess based on first P01)
    is_raining = rand < P01_t(1);
    
    for i = 1:n
        % Determine transition probability based on previous state
        if is_raining
            prob_wet = P11_t(i);
        else
            prob_wet = P01_t(i);
        end
        
        % Update State
        if rand < prob_wet
            is_raining = true;
            % Generate Amount (Gamma Distribution)
            % gamrnd(shape, scale)
            precip(i) = gamrnd(Alpha_t(i), Beta_t(i));
            
            % Enforce threshold floor (optional, improves realism for trace amounts)
            if precip(i) < coeffs.wet_threshold
                 precip(i) = coeffs.wet_threshold;
            end
        else
            is_raining = false;
            precip(i) = 0;
        end
    end
end