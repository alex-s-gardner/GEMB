function simulate_coeffs_disp(coeffs, struct_name)
% simulate_coeffs_disp Displays structure fields as executable MATLAB code.
%
%% Syntax
%
% simulate_coeffs_disp(coeffs, struct_name)
%
%% Description
%
% This utility function iterates through the fields of a given coefficient 
% structure and prints them to the command window. 
% The output is formatted as valid MATLAB assignment statements (e.g., 
% struct.field = value;), making it easy to copy, paste, and reuse parameter 
% sets for debugging or configuration files.
%
%% Inputs
%
%  coeffs       : struct       Structure containing coefficient values (scalars or vectors).
%  struct_name  : string       (Optional) Name of the structure variable to be used in the output string. Defaults to "coeffs".
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

if nargin < 2
    struct_name  = "coeffs";
end

keys = fieldnames(coeffs);
for i = 1:length(keys)
    key = keys{i};
    if length(coeffs.(key)) > 1 
        if size(coeffs.(key), 1) == 1
            disp(struct_name + "." + key + " = [" + sprintf('%0.4f ', coeffs.(key)) + "];")
        else
            disp(struct_name + "." + key + " = [" + sprintf('%0.4f ', coeffs.(key)) + "]';")
        end
    else
        disp(struct_name  + "." + key + " = " + sprintf('%0.4f ', coeffs.(key)) + ";")
    end
end
disp(" ")
end
