function dlw_delta = simulate_longwave_irradiance_delta(dec_year, coeff)
    n = length(dec_year);
    dlw_delta = zeros(n, 1);
    
    % 1. Component Selection (1, 2, or 3)
    % Cumulative probability used for efficient multi-component sampling
    cumP = cumsum(coeff.P);
    comp_selector = rand(n, 1);
    
    % 2. Sample and Truncate
    for i = 1:n
        % Determine which component to sample
        idx = find(comp_selector(i) <= cumP, 1, 'first');
        
        % Draw until value is within observed 'soft' bounds
        % This prevents the 'spread' issue while keeping the 3-peak shape
        val = coeff.mu(idx) + coeff.sigma(idx) * randn();
        while val < coeff.prctile_bounds(1) || val > coeff.prctile_bounds(2)
             val = coeff.mu(idx) + coeff.sigma(idx) * randn();
        end
        dlw_delta(i) = val;
    end
end