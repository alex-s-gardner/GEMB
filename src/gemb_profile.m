function Profile = gemb_profile(OutData, index)
% gemb_profile generates a tabulated Profile from a gemb output structure.
%
%% Syntax
% 
%  Profile = gemb_profile(OutData)
%  Profile = gemb_profile(OutData, index)
%
%% Description
%
% Profile = gemb_profile(OutData) creates a Profile table containing column
% properties from the last time step OutData from the gemb function. 
% 
% Profile = gemb_profile(OutData, index) specifies an integer timestep index 
% in the range of 1 to the length of OutData.time. 
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
    index (1,1) {mustBeInteger} =numel(OutData.time)
end

assert(index>0,'Input index must be greater than zero.')
assert(index<=numel(OutData.time),'Input index must be greater than zero.')

%%

dz               = OutData.dz(:,index); 
z_center         = dz2z(dz); 
temperature      = OutData.temperature(:, index); 
density          = OutData.density(:, index); 
water            = OutData.water(:, index); 
grain_radius     = OutData.grain_radius(:, index); 
grain_dendricity = OutData.grain_dendricity(:, index); 
grain_sphericity = OutData.grain_sphericity(:, index); 
albedo           = OutData.albedo(:, index); 
albedo_diffuse   = OutData.albedo_diffuse(:, index); 

Profile = table(z_center, dz, temperature, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse);

end

