function z_center = dz2z(dz)
% dz2z returns a center coordinates from GEMB column spacings. 
% 
%% Syntax 
% 
%  z_center = dz2z(dz)
%
%% Description
% 
% z_center = dz2z(dz) the center heights of a GEMB column
% each grid cell in meters.
% 
%% Example 1: Visualize the initial column: 
% Initialize a column and plot its grid spacing. A mean air temperature is
% defined below because it is used by model_initialize_column to set the
% initial temperature of the snow/firn/ice column. 
%
%   % Initialize parameters: 
%   ModelParam = model_initialize_parameters();
%   ClimateForcing.T_air_mean = 253.15; % -20 C
%  
%   % Initialize Column: 
%   [~, dz] = model_initialize_column(ModelParam, ClimateForcing);
%
%   % Get height column corresponding to dz: 
%   z_center = dz2z(dz);
%
%   % Plot 
%   plot(dz,z_center,'o-')
%   xlabel 'Vertical grid spacing (m)' 
%   ylabel 'Column height (m)' 
% 
%% Example 2: Visualize GEMB results: 
% Create a Hovmöller diagram of snow/firn/ice column temperature using
% pcolor. 
%
%   % Initialize model parameters:
%   ModelParam = model_initialize_parameters();
%   ModelParam.output_frequency = "daily"; 
%   
%   % Generate sample data: 
%   time_step_hours = 3;
%   ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);
%   
%   % Initialize grid:
%   [T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(ModelParam, ClimateForcing);
%   
%   % Run GEMB: 
%   OutData = gemb(T, dz, d, W, re, gdn, gsp, a, a_diffuse, ClimateForcing, ModelParam);
%   
%   % Get a 2D matrix of grid cell centers: 
%   z_center = dz2z(OutData.dz);
%   
%   % As of the writing of this example, OutData.time is all NaN, so we
%   will plot "time" as output timesteps and convert it to 2D so pcolor can
%   plot it: 
%   time_2D = repmat(1:numel(OutData.time),size(OutData.T,1),1);
%   
%   figure
%   pcolor(time_2D,z_center,OutData.T)
%   shading interp
%   clim([250 270])
%   ylabel 'Column height (m)'
%   xlabel 'Time step'
%   ylim([-10 1])
%   cb = colorbar;
%   ylabel(cb,'Temperature (K)')
%   cmocean thermal % optional colormap
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

arguments
    dz (:,:) 
end

% The bottom of the top grid cell is located at -dz(first_finite_value), so
% the center of the top grid cell is half a grid cell above that:
z_center = -cumsum(dz,'omitnan') + get_first_finite_data(dz)/2; 

end

function first_finite_data = get_first_finite_data(A)
% get_first_finite_data returns the top row of finite data in matrix A.
% This function is needed because OutData.dz is a 2D matrix that may have
% many rows of NaNs at the top

% Create a logical mask of finite values
mask = isfinite(A);

% Use max to find the index of the first 'true' (1) in each column
[has_finite, first_idx] = max(mask, [], 1);

% Filter out columns that have no finite values at all
first_idx(~has_finite) = NaN; 

% Extract the values using linear indexing
valid_cols = find(has_finite);
first_finite_data = nan(1, size(A, 2));
first_finite_data(valid_cols) = A(sub2ind(size(A), first_idx(valid_cols), valid_cols));

end