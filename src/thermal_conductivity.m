function K = thermal_conductivity(temperature, density, ModelParam)
% thermal_conductivity computes the thermal conductivity profile for snow, 
% firn, and ice based on density and temperature.
%
%% Syntax
%
% K = thermal_conductivity(temperature, density, ModelParam)
%
%% Description
%
% This function calculates the effective thermal conductivity of the model 
% column . It differentiates between snow/firn and solid glacier ice 
% based on a density threshold defined in ModelParam.density_ice.
%
% The calculation uses distinct parameterizations for each regime:
% 1. Snow and Firn (density < density_ice): Thermal conductivity is calculated 
%    as a function of density using empirical regressions. The user can select 
%    between the formulations of:
%    * Sturm et al. (1997): K = 0.138 - 1.01e-3*density + 3.233e-6*density^2.
%    * Calonne et al. (2011): K = 0.024 - 1.23e-4*density + 2.5e-6*density^2.
% 2. Glacier Ice (density >= density_ice): Thermal conductivity is dominated 
%    by phonon transport and is calculated as a function of temperature:
%    K = 9.828 * exp(-5.7e-3 * temperature).
%
%% Inputs
%
%  temperature                      : K            Grid cell temperature (vector).
%  density                          : kg m^-3      Grid cell density (vector).
%  ModelParam                       : struct       Model parameter structure:
%    .density_ice                   : kg m^-3      Density threshold defining glacier ice.
%    .thermal_conductivity_method   : string       Parameterization choice: "Sturm" or "Calonne".
%
%% Outputs
%
%  K                                : W m^-1 K^-1  Vector of thermal conductivities.
%
%% References
% If you use GEMB, please cite the following:
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.
%
% Underlying physical parameterizations:
% Sturm, M., et al. (1997). The thermal conductivity of seasonal snow.
%   Journal of Glaciology.
% Calonne, N., et al. (2011). Thermal conductivity of snow...
%   Geophysical Research Letters.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

%% CONSTANTS & INITIALIZATION

% Tolerance to prevent floating point errors at the density threshold
d_tolerance  = 1e-11;

% Get number of grid cells
m = length(density);

% Initialize conductivity vector with zeros
K = zeros(m,1);

%% IDENTIFY SNOW VS ICE

% Create logical mask: True for snow/firn, False for ice
sfIdx = density < ModelParam.density_ice - d_tolerance ;

%% CALCULATE CONDUCTIVITY FOR SNOW/FIRN

% Use empirical density-based regressions
switch ModelParam.thermal_conductivity_method
    case "Calonne"
        % Parameterization from Calonne et al. (2011)
        % Often used for a wider range of snow microstructures
        K(sfIdx) = 0.024 - 1.23E-4 * density(sfIdx) + 2.5e-6 * (density(sfIdx).^2);
        
    case "Sturm"
        % Parameterization from Sturm et al. (1997) [Default]
        % Standard regression for seasonal snow
        K(sfIdx) = 0.138 - 1.01E-3 * density(sfIdx) + 3.233E-6 * (density(sfIdx).^2);
end

%% CALCULATE CONDUCTIVITY FOR ICE
% For densities >= ModelParam.density_ice, conductivity is dominated by temperature dependence.
% Formula typically attributed to Weller & Schwerdtfeger (1977) or similar
% standard glaciological relations.
% Note: ~sfIdx selects the inverse of the snow index (i.e., the ice cells)

K(~sfIdx) = 9.828 * exp(-5.7E-3 * temperature(~sfIdx));

end