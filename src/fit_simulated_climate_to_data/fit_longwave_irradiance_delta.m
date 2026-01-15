function coeff = fit_longwave_irradiance_delta(delta)
% fit_longwave_irradiance_delta 
%
%% Syntax
% 
%  coeff = fit_longwave_irradiance_delta(delta)
% 
%% Description
%
% coeff = fit_longwave_irradiance_delta(delta)
%
%% Example 
% 
%
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

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