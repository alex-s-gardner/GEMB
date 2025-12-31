function datenum_out = decyear2datenum(decyear)
    % DECYEAR2DATENUM Converts decimal year to MATLAB serial date number.
    %
    % Usage: datenum_out = decyear2datenum(2023.5)
    
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