function K = thermal_conductivity(d, T, dIce, thermal_conductivity_method)
% thermal_conductivity computes the thermal conductivity profile for snow, 
% firn, and ice based on density and temperature.
%
%% Syntax
%   K = thermal_conductivity(d, T, dIce, thermal_conductivity_method)
%
%% Description
%   Calculates thermal conductivity [W m-1 K-1] differentiating between 
%   snow/firn (density < dIce) and glacier ice (density >= dIce).
%   
%   For snow/firn, it uses empirical relationships based on density.
%   For ice, it uses a temperature-dependent relationship.
%
%% Inputs
%   d     : vector of grid cell densities [kg m-3]
%   T     : vector of grid cell temperatures [K]
%   dIce  : density threshold defining glacier ice (e.g., 910 or 917) [kg m-3]
%   thermal_conductivity_method : integer flag for snow conductivity parameterization:
%           1 = Sturm et al. (1997) [Default]
%           2 = Calonne et al. (2011)
%
%% Outputs
%   K     : vector of thermal conductivities [W m-1 K-1]
%
%% References
%   Sturm, M., et al. (1997). The thermal conductivity of seasonal snow.
%       Journal of Glaciology.
%   Calonne, N., et al. (2011). Thermal conductivity of snow...
%       Geophysical Research Letters.

    %% CONSTANTS & INITIALIZATION
    % Tolerance to prevent floating point errors at the density threshold
    Dtol = 1e-11;
    
    % Get number of grid cells
    m = length(d);
    
    % Initialize conductivity vector with zeros
    K = zeros(m,1);
    
    %% IDENTIFY SNOW VS ICE
    % Create logical mask: True for snow/firn, False for ice
    sfIdx = d < dIce - Dtol;
    
    %% CALCULATE CONDUCTIVITY FOR SNOW/FIRN
    % Use empirical density-based regressions
    if thermal_conductivity_method == 2
        % Parameterization from Calonne et al. (2011)
        % Often used for a wider range of snow microstructures
        K(sfIdx) = 0.024 - 1.23E-4 * d(sfIdx) + 2.5e-6 * (d(sfIdx).^2);
    else 
        % Parameterization from Sturm et al. (1997) [Default]
        % Standard regression for seasonal snow
        K(sfIdx) = 0.138 - 1.01E-3 * d(sfIdx) + 3.233E-6 * (d(sfIdx).^2);
    end
    
    %% CALCULATE CONDUCTIVITY FOR ICE
    % For densities >= dIce, conductivity is dominated by temperature dependence.
    % Formula typically attributed to Weller & Schwerdtfeger (1977) or similar
    % standard glaciological relations.
    % Note: ~sfIdx selects the inverse of the snow index (i.e., the ice cells)
    K(~sfIdx) = 9.828 * exp(-5.7E-3 * T(~sfIdx));
end