function S = combineStrucData_GEMB(S,LP,runIdx)
% combine input data into structure for passing to parfor loop
% 
%% Syntax 
% 
% 
%
%% Description
% 
% 
% 
%% Inputs
% 
% 
% 
%% Outputs
% 
% 
%% Documentation
% 
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
% 
%% References 
% If you use GEMB, please cite the following: 
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass 
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci. 
% Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

% transfer local parameters to S structure
fn = fields(LP);
for i = 1:length(fn)
    S.(fn{i}) = LP.(fn{i});
end

S.run_id = [S.run_prefix '_' sprintf('%06d', runIdx)];
