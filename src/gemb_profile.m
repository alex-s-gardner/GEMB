function Profile = gemb_profile(OutData, time_extract)
% gemb_profile generates a tabulated Profile from a gemb output structure. The Profile generated 
% by this function can be used to restart gemb from a saved state such as after a spinup run. 
%
%% Syntax
% 
%  Profile = gemb_profile(OutData)
%  Profile = gemb_profile(OutData, time_extract)
%
%% Description
%
% Profile = gemb_profile(OutData) creates a Profile table containing column
% properties from the last time step OutData from the gemb function. 
% 
% Profile = gemb_profile(OutData, time_extract) specifies a query datetime time_extract
% corresponding to the desired output profile. If time_extract is not specified,
% the last time step OutData is used. If time_extract does not exactly match any
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
    time_extract (1,1) datetime = OutData.time(end)
end

assert(time_extract>=OutData.time(1)  ,'time_extract cannot be before the first time step of OutData.')
assert(time_extract<=OutData.time(end),'time_extract cannot be after the last time step of OutData.')

%% Construct the table

% Get the index of the closest time step: 
if time_extract == OutData.time(end)
    % interp fails when length(OutData.time) == 1
    index = length(OutData.time);
else
    index = interp1(OutData.time, 1:numel(OutData.time), time_extract, "nearest"); 
end

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