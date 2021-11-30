function [dateN, Ta0, V0, dsw0, dlw0, P0, eAir0, pAir0, dt, LP] = ...
    loadMetData(cldFrac, site, dataFreq, saveFreq)
%% LOAD IN PR-PROCESSED MODEL FORCING

% Description:
% loads in model forcing from either a series of AWS measurments collected 
% over the Devon Island Ice Cap or creates synthetic data using 
% 'metDataCreate.m'

%% Site specifications
% file name of .mat file containing meterological data saved in the
% specified structure format

% fileName = 'D_1_1hr_met_data';
% TITLE = 'Site D1 2006' %Location title to be included in figures
% 
% site location
%   LP.lat = 75.3399;      % decimal deg North
%   LP.lon = 277.3237;     % decimal deg East
%   LP.elev = 1802;        % m above mean sea level

% REMOVE SWITCH @ LATER DATE
% switch made for sites D1 (1), D2(2) and D3(3)

% * dateN : date/time [UTC]
% * Ta: 2m air temperature [K]
% * V: wind speed [m s^-^1]
% * dsw: downward shortwave radiation flux [W m^-^2]
% * dlw: downward longwave radiation flux [W m^-^2] 
% * P: precipitation [kg m^-^2]
% * eAir: screen level vapor pressure [Pa]
% * eAir: screen level vapor pressure [Pa]
% * dt: time step of met data [s]
% * LP.Tmean: mean annual temperature [K]
% * LP.C = mean annual snow accumulation [kg m-2 yr-1]
% * LP.elev: surface elevation [m a.s.l.]
% * LP.Tz: height above ground at which temperature (T) was sampled [m]
% * LP.Vz: height above ground at which wind (V) eas sampled [m]
% * LP.elev: surface elevation [m a.s.l.]

% set local parameteres
% specify mean annual temperature [K]
LP.Tmean = -15 + 273.15;

% height above ground at which temperature (T) and wind (v) were measured [m]
LP.Vz = 2;
LP.Tz = 2;

% specify mean annual snow accumulation [kg m-2 yr-1]
LP.C = 150;


switch site
    case 0
        LP.lat = 75.3399;    
        LP.lon = 277.3237;    
        LP.elev = 1802;
        addpath(fullfile('..','Input','Artificial'))
        met = metDataCreate(LP.lat, LP.lon, LP.elev);
                
    case 1     
        if dataFreq == 1;
            fileName = 'D_1_daily_met_data';
        else
            fileName = 'D_1_1hr_met_data';
        end
        
        LP.lat = 75.3399;    
        LP.lon = 277.3237;    
        LP.elev = 1802;
 
    case 2
        if dataFreq == 1;
            fileName = 'D_2_daily_met_data';
        else
            fileName = 'D_2_1hr_met_data';
        end

        LP.lat = 75.1789;     
        LP.lon = 277.2248;    
        LP.elev = 1415;
        
    case 3 
        if dataFreq == 1;
            fileName = 'D_3_daily_met_data';
        else
            fileName = 'D_3_1hr_met_data';
        end

        LP.lat = 75.0139;     
        LP.lon = 277.1236;    
        LP.elev = 994; 
end

if site ~= 0                    % don't load met data if auto generated
    load (fullfile('..','Input',fileName))
end

dateN = met.date;               % [UTC]
Ta0 = met.T_air + 273.15;       % [K]
V0 = met.V;                     % [m s-1]
dsw0 = met.dswrf;             % [w m-2]
RH0 = met.RH;                   % [%]
P0 = met.Precip;                % [mm w.e.] == [kg m-2]

% if data is provided hourly and not artifically generated
if dataFreq == 2 && site ~= 0                      
    model.z_measure  = ...              % snow depth [m] measured with UDG
        met.z_snow(saveFreq:saveFreq:end)';      
    model.Q_net_measure = ...           % measured net radiation [W m-2]
        met.Rad_net(saveFreq:saveFreq:end)';  
end

% saturated vapor pressure according to [Pa]
eSat = 100 .* 10.^(-7.90298 .* (373.16 ./ Ta0 - 1)... 
    + 5.02808 .* log10(373.16 ./ Ta0)...
    - 1.3816E-7 .* (10.^(11.344 .* (1 - Ta0 ./ 373.16))-1)...
    + 8.1328E-3*(10.^(-3.49149.*(373.16./Ta0-1))-1)...
    + log10(1013.246)); 

% calculate screen level vapor pressure [Pa]
eAir0 = eSat .* (RH0/100);

% calculate incoming longwave radiation [W m-2]
dlw0 = longwaveDn(Ta0, eAir0, cldFrac);

% determine data frequency
dt = (dateN(2)-dateN(1)) * 86400; % data frequency [s]
dt = round(dt*100) /100;                     % round dt to milliseconds


% calculated air pressure using the Barometric Formula [Pa]
% M_air = 0.029; % molecular mass of air [kg/mol]
% R = 8.314; % Universal gas constant for air [N·m /(mol·K)]
Ta_sea = Ta0 + (LP.elev * 0.0065);  % temperature at sea level

pAir0 = 101325 * exp((-0.029 * 9.81 * LP.elev) ./ (8.314 * Ta_sea));

%% Check met data quality
% make sure the the data frequency devides evenly into a day
if rem(86400, dt) > 0
    error('data frequency does devide evenly into a day')
end

% check for missing data
if sum(isnan(dateN)) > 0
    warning('some or all date numbers missing')
elseif sum(isnan(Ta0)) > 0
    error('some or all air temperature data missing')
elseif sum(isnan(V0)) > 0
    error('some or all wind speed data missing')
elseif sum(isnan(dsw0)) > 0
    error('some or all incoming shortwave data missing')
elseif sum(isnan(P0)) > 0
    error('some or all precipitation data missing')
end

% check for erroneous values
if sum(Ta0 > 333.15) > 0 || sum(Ta0 < 208.15)
    warning(['air T > 60 or <-65 °C)'...
        'temperature must be specified in °C'])
elseif sum(V0 < 0) > 0 || sum(V0 > 45) > 0
    warning('wind speed in < 0 or < 45 m/s')
elseif sum(dsw0 < 0) > 0 || sum(dsw0 > 1400) > 0
    warning('solar radiation is < 0 or > 1400')
elseif sum(RH0 < 0) > 0 || sum(RH0 > 100) > 0
    warning('RH is < 0 or > 100')
elseif sum(P0 < 0) > 0 || sum(P0 > 100) > 0
    warning('Pricipitation is < 0 or > 100 mm w.e.')
end

end

function dlw = longwaveDn(Ta, eAir, cldFrac)

% Description: calculates clear-sky incoming longwave radiation [W m-2] 
% based on screen level temperature (Ta)[K] and screen level vapor
% pressure [eAir]
% 
% equations are based on 2 years of on-ice measurements taken over the
% Greenland Ice Cap
 
% Reference:
% Konzelmann, T., R. S. W. Vandewal, W. Greuell, R. Bintanja, E. A. C.
% Henneken, and A. Abeouchi, 1994: Parameterization of Global and Longwave
% Incoming Radiation for the Greenland Ice-Sheet. Global and Planetary
% Change, 9, 143-164.
 
%% clear sky longwave emittance 
% m = 8;
% % clear sky emittance 
% e_cs = 0.23 + b .* (eAir ./ Ta).^(1/m); 
% % downward longwave radiation [w m-2]
% dlw = e_cs .* (5.67E-8 * Ta.^4);

%% daily average longwave emittance including clouds
dlw = ((0.23 + 0.483.*(eAir./Ta).^(1/8)).*(1-cldFrac^3)+...
    0.963*cldFrac^3).* 5.67E-8 .* Ta.^4;
end