function A_surface = surface_timeseries(A)
% surface_timeseries returns the top row of finite data in a matrix.
% 
%% Syntax
%
%  A_surface = surface_timeseries(A)
%
%% Description 
%
% A_surface = surface_timeseries(A) for input matrix A of dimensions MxN
% returns a 1xN vector A_surface of the uppermost finite values of the 
% matrix A.
% 
%% Example 
% 
%   % Generate sample data: 
%   time_step_hours = 3;
%   ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);
%   
%   % Initialize model parameters:
%   ModelParam = model_initialize_parameters(output_frequency="daily");
%   
%   % Initialize grid:
%   [temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse] = model_initialize_column(ModelParam, ClimateForcing);
%  
%   % Run GEMB: 
%   OutData = gemb(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, ClimateForcing, ModelParam);
%   
%   % Get a time series of skin temperature: 
%   temperature_skin = surface_timeseries(OutData.temperature); 
% 
%   % Create a datetime array to make it easy: 
%   dates_datetime = datetime(OutData.dates,'convertfrom','datenum'); 
% 
%   % Plot the timeseries: 
%   plot(dates_datetime, temperature_skin)
%   ylabel 'Skin temperature (K)'
% 
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
% 
% See also dz2z. 


% Input checks: 
assert(ismatrix(A) & ~isscalar(A) & size(A,1)>1,'Input matrix A must be a 2D matrix.')

% Create a logical mask of finite values
mask = isfinite(A);

% Use max to find the index of the first 'true' (1) in each column
[has_finite, first_idx] = max(mask, [], 1);

% Filter out columns that have no finite values at all
first_idx(~has_finite) = NaN; 

% Extract the values using linear indexing
valid_cols = find(has_finite);
A_surface = nan(1, size(A, 2));
A_surface(valid_cols) = A(sub2ind(size(A), first_idx(valid_cols), valid_cols));

end