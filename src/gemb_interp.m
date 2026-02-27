function A_regularized = gemb_interp(Z_center,A,Profile,interp_method)
% gemb_interp regularizes the MxN gemb output onto a consistent vertical grid. 
% This function is necessary because the vertical spacing of GEMB output
% evolves with every timestep. 
%
%% Syntax
%
%  A_regularized = gemb_interp(Z_center,A,Profile)
%  A_regularized = gemb_interp(Z_center,A,Profile,interp_method)
%
%% Description 
%
% A_regularized = gemb_interp(Z_center,A,Profile) regrids the MxN (depth x
% time) data A onto an M_regularized x N grid, where the vertical postings
% of the new grid correspond to the z_center field of the table Profile.
% Interpolation is performed vertically for each timestep. 
%
% A_regularized = gemb_interp(Z_center,A,Profile,interp_method) specifies
% a 1D interpolation method. Default is "linear". 
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

%% Parse inuts 

arguments
    Z_center (:,:)  {mustBeNumeric}
    A        (:,:)  {mustBeNumeric}
    Profile  (:,10) table
    interp_method {mustBeText} = "linear"
end

assert(isequal(size(Z_center),size(A)),"Dimensions of Z_center and A must agree.")
assert(size(Z_center,1)>1, "Inputs Z_center and A must contain multiple rows representing profile depth.")

%% Interpolate

isf = isfinite(A); 

% Preallocate output: 
A_regularized = NaN(height(Profile),size(A,2));

% Loop through each timestep:
for k = 1:size(A,2)

    A_regularized(:,k) = interp1(Z_center(isf(:,k),k),A(isf(:,k),k),...
        Profile.z_center,interp_method,"extrap"); 

end

end