function coeffs = fit_precipitation(dec_year, precip)
% FIT_PRECIPITATION Fits a seasonal Markov-Gamma model to hourly precipitation.
%
%   INPUTS:
%       dec_year - N x 1 vector of decimal years (e.g., 2021.45)
%       precip   - N x 1 vector of precipitation amounts (mm)
%
%   OUTPUT:
%       coeffs   - Struct containing harmonic coefficients for the model
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

    % 1. Constants and Setup
    wet_threshold = 0.1; % Minimum mm to count as "wet"
    
    % Extract fractional year (0 to 1) for seasonality
    t_season = mod(dec_year, 1);
    
    % Determine Wet/Dry states
    is_wet = precip >= wet_threshold;
    
    % 2. Calculate Monthly Statistics (12 bins)
    % We bin by month to get robust estimates, then fit harmonics to the bins.
    months = ceil(t_season * 12);
    months(months == 0) = 1; % Handle edge case exactly at 0.0
    
    stats = zeros(12, 4); % Columns: [P01, P11, Alpha, Beta]
    
    for m = 1:12
        idx = (months == m);
        data_m = is_wet(idx);
        precip_m = precip(idx);
        
        if sum(idx) < 2
            continue; 
        end
        
        % --- Markov Transition Probabilities ---
        % Current state vs Next state (shift by 1)
        % Note: This simple binning ignores the transition between months, 
        % which is acceptable for parameter estimation.
        current = data_m(1:end-1);
        next    = data_m(2:end);
        
        % P01: Probability Dry -> Wet (Count 0->1 / Count 0)
        idx_0 = (current == 0);
        if sum(idx_0) > 0
            stats(m, 1) = sum(next(idx_0) == 1) / sum(idx_0);
        end
        
        % P11: Probability Wet -> Wet (Count 1->1 / Count 1)
        idx_1 = (current == 1);
        if sum(idx_1) > 0
            stats(m, 2) = sum(next(idx_1) == 1) / sum(idx_1);
        end
        
        % --- Gamma Distribution Parameters (Method of Moments) ---
        wet_amts = precip_m(precip_m >= wet_threshold);
        if length(wet_amts) > 5
            mu = mean(wet_amts);
            v  = var(wet_amts);
            
            % Gamma params: Alpha (Shape) = mean^2 / var, Beta (Scale) = var / mean
            if v > 0
                stats(m, 3) = mu^2 / v; % Alpha
                stats(m, 4) = v / mu;   % Beta
            end
        end
    end
    
    % 3. Fit Harmonic Curves
    % Model: Y = C1 + C2*sin(2*pi*t) + C3*cos(2*pi*t)
    % We fit this to the 12 monthly values.
    
    t_bins = (0.5:1:11.5)' / 12; % Midpoints of months
    X = [ones(12,1), sin(2*pi*t_bins), cos(2*pi*t_bins)];
    
    % Solve least squares: B = (X'X)^-1 X'Y
    % We do this for all 4 parameters at once
    B = X \ stats; 
    
    % Pack into structure
    coeffs.P01_harmonics = B(:, 1);
    coeffs.P11_harmonics = B(:, 2);
    coeffs.Alpha_harmonics = B(:, 3);
    coeffs.Beta_harmonics  = B(:, 4);
    
    % Store threshold for simulation
    coeffs.wet_threshold = wet_threshold;
end