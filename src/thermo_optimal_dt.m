function dt = thermo_optimal_dt(dz, d, CI, K, global_dt_or_dt_divisors)
    % determine minimum acceptable delta t (diffusion number > 1/2) [s]
    % 1. Calculate the theoretical limit for every single grid cell
    %    (Using 0.5 is the absolute limit; we usually aim lower for safety)
    stability_limit_per_cell = 0.5 * (d .* CI .* dz.^2) ./ K;
    
    % 2. Find the bottleneck: The smallest allowable dt in the entire column
    max_safe_dt = min(stability_limit_per_cell);
    
    % 3. Apply a Safety Factor (0.9 or 0.8 is standard)
    %    This accounts for floating point errors or slight non-linearities in K
    dt_target = max_safe_dt * 0.8;
    
    % 4. (Optional) Sanity check to prevent extremely small steps if bad data enters
    if dt_target < 1e-4
        warning('Timestep is extremely small (%e). Check for near-zero dz layers.', dt_target);
    end
    
    % 5. Fit this target into your input data frequency (global_dt)
    %    Find the largest divisor of global_dt that is <= dt_target  
    if isscalar(global_dt_or_dt_divisors)
        dt_divisors = fast_divisors(global_dt_or_dt_divisors * 10000)/10000; % ClimateForcingStep.dt is in seconds fastDivisors is a subfunction below. 
    else
        dt_divisors = global_dt_or_dt_divisors;
    end

    dt = dt_divisors(find(dt_divisors <= dt_target, 1, 'last'));
    
    if isempty(dt)
        warning("thermo dt_target < all dt_divisors, setting thermo == to the smallest dt_divisors... this may make termo diffusion unstable")
        dt = dt_divisors(1); % Fallback to smallest possible step
    end

end