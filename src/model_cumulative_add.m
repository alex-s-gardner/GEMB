function OutCum = model_cumulative_add(M, R, F, EC, Ra, M_added, ...
    sw_net, lw_net, shf, lhf, ulw, compaction_dens, compaction_melt, ...
    d, a, re, dz, ModelParam, OutCum)
% model_cumulative_add updates cumulative variables for model output.
%
%% Syntax
%
% OutCum = model_cumulative_add(OutCum, M, R, F, EC, Ra, M_added, ...
%    sw_net, lw_net, shf, lhf, ulw, compaction_dens, compaction_melt, ...
%    d, a, re, dz, ModelParam)
%
%% Description
%
% This function updates the tracking structure `OutCum` by adding the current
% time-step's values to the running totals. It replaces the use of `eval()`
% in earlier versions of GEMB to significantly improve runtime performance.
%
% It performs two main tasks:
% 1. Calculates derived variables for the current state (e.g., Q_net, FAC, 
%    surface properties d1, a1, re1).
% 2. Explicitly sums these values into the `OutCum` structure fields.
%
%% Inputs
%
%  OutCum           : struct       Structure containing cumulative variables from previous steps.
%  M                : kg m^-2      Melt mass.
%  R                : kg m^-2      Runoff mass.
%  F                : kg m^-2      Refrozen mass.
%  EC               : kg m^-2      Evaporation/Condensation mass.
%  Ra               : kg m^-2      Rain mass.
%  M_added          : kg m^-2      Mass added/removed by layer management.
%  sw_net           : W m^-2       Net shortwave radiation.
%  lw_net           : W m^-2       Net longwave radiation.
%  shf              : W m^-2       Sensible heat flux.
%  lhf              : W m^-2       Latent heat flux.
%  ulw              : W m^-2       Upward longwave radiation.
%  compaction_dens  : m            Compaction due to densification.
%  compaction_melt  : m            Compaction due to melt.
%  d                : kg m^-3      Density profile.
%  a                : fraction     Albedo profile.
%  re               : mm           Grain radius profile.
%  dz               : m            Layer thickness profile.
%  ModelParam       : struct       Model parameters (needs .density_ice).
%
%% Outputs
%
%  OutCum           : struct       Updated cumulative structure with incremented .count.
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% 1. Calculate derived variables for output
d1    = d(1);
a1    = a(1);
re1   = re(1);
Q_net = sw_net + lw_net + shf + lhf;

% Firn Air Content (FAC) [m]
% Defined as the integrated column thickness of air equivalent.
% FAC = sum(dz * (rho_ice - rho) / rho_ice) for rho < rho_ice
% Note: The original implementation divided by 1000 instead of ModelParam.density_ice?
% Preserving original logic: sum(dz.*(ModelParam.density_ice - min(d,ModelParam.density_ice)))/1000;
FAC = sum(dz .* (ModelParam.density_ice - min(d, ModelParam.density_ice))) / 1000;

% 2. Explicitly accumulate values
% Using explicit assignment is significantly faster than dynamic field access
OutCum.R               = OutCum.R + R;
OutCum.M               = OutCum.M + M;
OutCum.F               = OutCum.F + F;
OutCum.EC              = OutCum.EC + EC;
OutCum.Ra              = OutCum.Ra + Ra;
OutCum.M_added         = OutCum.M_added + M_added;

OutCum.sw_net          = OutCum.sw_net + sw_net;
OutCum.lw_net          = OutCum.lw_net + lw_net;
OutCum.shf             = OutCum.shf + shf;
OutCum.lhf             = OutCum.lhf + lhf;
OutCum.ulw             = OutCum.ulw + ulw;
OutCum.Q_net           = OutCum.Q_net + Q_net;

OutCum.a1              = OutCum.a1 + a1;
OutCum.re1             = OutCum.re1 + re1;
OutCum.d1              = OutCum.d1 + d1;

OutCum.compaction_dens = OutCum.compaction_dens + compaction_dens;
OutCum.compaction_melt = OutCum.compaction_melt + compaction_melt;
OutCum.FAC             = OutCum.FAC + FAC;

% Increment the counter
OutCum.count = OutCum.count + 1;

end