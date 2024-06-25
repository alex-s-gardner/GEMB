function jd = juliandate(varargin)
% This sub function is provided in case juliandate does not come with your 
% distribution of Matlab
[year month day hour min sec] = datevec(datenum(varargin{:}));
idx = month <= 2;
year(idx) = year(idx)-1;
month(idx) = month(idx)+12;
jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
    floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
    (hour + min/60 + sec/3600)/24;
end
