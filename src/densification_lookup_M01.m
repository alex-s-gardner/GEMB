function M01 = densification_lookup_M01(densification_coeffs_M01)
% densification_lookup_M01 returns calibrated coefficients of a densification model. 
% 
%% Model Description
% See Section 4.2 Firn air content of Gardner et al (2023). The model is based
% on Ligtenberg et al (2011)'s Equations 8 and 9, where the output of this
% function matches the coefficients MO_550 and MO_830.
%
%  MO_density = MO_density_offset - MO_density_slope * log(accumulation_rate)
%
%% Syntax
% 
%  M01 = densification_lookup_M01(densification_coeffs_M01) 
%
%% Description
% 
% M01 = densification_lookup_M01(densification_coeffs_M01) returns
% coefficients in the form: 
% 
%   M01 = [M0_550_offset M0_550_slope M0_830_offset M0_830_slope] 
%
% or 
% 
%   M01 = [M0_550_offset M0_550_slope M0_830_offset M0_830_slope;
%          M1_550_offset M1_550_slope M1_830_offset M1_830_slope] 
% 
%% Example 
% 
%   densification_lookup_M01("Gre_ERA5_GS_SW0")
%   ans =
%       1.3566    0.1350    1.8705    0.2290
%       1.4318    0.1055    2.0453    0.2137
%
%% References 
%
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
%
% Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R.: An improved 
% semi-empirical model for the densification of Antarctic firn, The Cryosphere, 
% 5, 809–819, https://doi.org/10.5194/tc-5-809-2011, 2011. 
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

switch densification_coeffs_M01
    % ------------------------ Antarctic -------------------------
    % ERA5 new albedo_method="GardnerSharp", shortwave_absorption_method=0
    case "Ant_ERA5_GS_SW0"
        M01 = [[1.5131, 0.1317, 0.1317, 0.2158]; [1.8422, 0.1688, 2.4979, 0.3225]];

    % ERA5 v4 (Paolo et al., 2023)
    case "Ant_ERA5v4_Paolo23"    
        M01 = [2.84, 0.32, 3.10, 0.37];

    % ERA5 new albedo_method="BruneLeFebre", shortwave_absorption_method=1
    case "Ant_ERA5_BF_SW1"
        M01 = [2.2191, 0.2301, 2.2917, 0.2710]; 

    % RACMO callibration, default (Gardner et al., 2023)
    case "Ant_RACMO_GS_SW0"
        M01 = [1.6383, 0.1691, 1.9991, 0.2414];

    %  Ligtenberg and others (2011), Antarctica
    case "Ant_Ligtenberg"
        M01 = [1.435, 0.151, 2.366, 0.293];

    % ------------------------- Greenland ------------------------
    % ERA5 new albedo_method="GardnerSharp", shortwave_absorption_method=0, bare ice
    case "Gre_ERA5_GS_SW0"
        M01 = [[1.3566, 0.1350, 1.8705, 0.2290]; [1.4318, 0.1055, 2.0453, 0.2137]];

    % RACMO callibration, default (Gardner et al., 2023)
    case "Gre_RACMO_GS_SW0"
        M01 = [1.2691, 0.1184, 1.9983, 0.2511];

    %  ismember(albedo_method,["GreuellKonzelmann","Bougamont2005"]) && shortwave_absorption_method>0
    case "Gre_RACMO_GB_SW1"
        M01 = [1.7834, 0.1409, 1.9260, 0.1527];
   
    % Kuipers Munneke and others (2015) [semi-empirical], Greenland
    case "Gre_KuipersMunneke"
        M01 = [1.042, 0.0916, 1.734, 0.2039];

    otherwise
        error("Unrecognized densification coefficients.")
end

