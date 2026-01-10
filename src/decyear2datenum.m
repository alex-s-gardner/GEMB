function datenum_out = decyear2datenum(decyear)
% DECYEAR2DATENUM Converts decimal year to MATLAB serial date number.
%
% Usage: datenum_out = decyear2datenum(2023.5)
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 

% 1. Extract the integer year
year_part = floor(decyear);

% 2. Calculate the start of the current year and the start of the next year
% This automatically accounts for leap years.
start_of_year = datenum(year_part, 1, 1);
start_of_next_year = datenum(year_part + 1, 1, 1);

% 3. Determine the total number of days in this specific year
days_in_year = start_of_next_year - start_of_year;

% 4. Calculate the fractional part of the year
fractional_year = decyear - year_part;

% 5. Compute the final date number
datenum_out = start_of_year + (fractional_year .* days_in_year);
end