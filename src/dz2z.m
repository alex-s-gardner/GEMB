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
%   ModelParam = model_initialize_parameters;
%   ClimateForcing.temperature_air_mean = 253.15; % -20 C
%  
%   % Initialize Column: 
%   Profile = model_initialize_column(ModelParam, ClimateForcing);
%
%   % Get height column corresponding to dz: 
%   z_center = dz2z(Profile.dz);
%
%   % Plot 
%   plot(Profile.dz,z_center,'o-')
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
%   Profile = model_initialize_column(ModelParam, ClimateForcing);
%   
%   % Run GEMB: 
%   OutData = gemb(Profile, ClimateForcing, ModelParam);
%   
%   % Get a 2D matrix of grid cell centers: 
%   z_center = dz2z(OutData.dz);
%   
%   % Convert time to 2D so pcolor can plot it: 
%   time_2D = repmat(OutData.dates,size(OutData.temperature,1),1);
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
% 
% See also surface_timeseries.

arguments
    dz (:,:) 
end

% The bottom of the top grid cell is located at -dz(first_finite_value), so
% the center of the top grid cell is half a grid cell above that:
z_center = -cumsum(dz,'omitnan') + surface_timeseries(dz)/2; 

end
