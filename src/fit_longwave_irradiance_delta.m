function coeff = fit_longwave_irradiance_delta(delta)
    % Clean data
    delta = delta(:);
    delta = delta(~isnan(delta));
    
    % Fit a 3-component Gaussian Mixture Model (Trimodal)
    options = statset('MaxIter', 1500, 'TolFun', 1e-6);
    % We use 3 components to better define the 'shoulder' and 'tails'
    gmm = fitgmdist(delta, 2, 'Options', options, 'RegularizationValue', 0.01);
    
    % Store the parameters
    coeff.mu    = gmm.mu;             
    coeff.sigma = sqrt(squeeze(gmm.Sigma)); 
    coeff.P     = gmm.ComponentProportion;  
    
    % Statistical range anchors to prevent range inflation
    coeff.limits = [min(delta), max(delta)];
    coeff.prctile_bounds = prctile(delta, [0.5, 95]); % Soft limits
end