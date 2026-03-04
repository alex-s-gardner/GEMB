function Profile = gemb_profile(OutData, time_i)
% gemb_profile generates a tabulated Profile from a gemb output structure.
%
%% Syntax
% 
%  Profile = gemb_profile(OutData)
%  Profile = gemb_profile(OutData, time_i)
%
%% Description
%
% Profile = gemb_profile(OutData) creates a Profile table containing column
% properties from the last time step OutData from the gemb function. 
% 
% Profile = gemb_profile(OutData, time_i) specifies a query datetime time_i
% corresponding to the desired output profile. If time_i is not specified,
% the last time step OutData is used. If time_i does not exactly match any
% elements in OutData.time, the profile corresponding to the nearest time step
% will be returned. 
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
% See also model_initialize_profile. 

%% Check inputs: 

arguments
    OutData struct 
    time_i (1,1) datetime =OutData.time(end)
end

assert(time_i>=OutData.time(1)  ,'time_i cannot be before the first time step of OutData.')
assert(time_i<=OutData.time(end),'time_i cannot be after the last time step of OutData.')

%% Construct the table

% Get the index of the closest time step: 
index = interp1(OutData.time, 1:numel(OutData.time), time_i, "nearest"); 

isf              = isfinite(OutData.dz(:,index)); 
dz               = OutData.dz(isf,index); 
z_center         = dz2z(dz); 
temperature      = OutData.temperature(isf, index); 
density          = OutData.density(isf, index); 
water            = OutData.water(isf, index); 
grain_radius     = OutData.grain_radius(isf, index); 
grain_dendricity = OutData.grain_dendricity(isf, index); 
grain_sphericity = OutData.grain_sphericity(isf, index); 
albedo           = OutData.albedo(isf, index); 
albedo_diffuse   = OutData.albedo_diffuse(isf, index); 

Profile = table(z_center, dz, temperature, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse);

end

