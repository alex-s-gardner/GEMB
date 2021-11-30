function [T, dz, d, W, a, re, gdn, gsp] = ...
    accumulation (T_air, T, dz, d, P, W, dz_min, a, a_SNOW, re, gdn, gsp)

%% Adds precipitation and deposition to the model grid

% Author: Alex Gardner, University of Alberta
% Date last modified: JAN, 2008

% Description:
% adjusts the properties of the top grid cell to account for accumulation
% T_air & T = Air and top grid cell temperatures [K]
% dz = topgrid cell length [m]
% d = density of top grid gell [kg m-3]
% P = precipitation [mm w.e.] or [kg m-3]
% re = effective grain radius [mm]
% gdn = grain dentricity
% gsp = grain sphericity

%% MAIN FUNCTION
% specify constants
dIce = 910;     % density of ice [kg m-3]
dSnow = 150;    % density of snow [kg m-3]
reNew = 0.1;    % new snow grain size [mm]
gdnNew = 1;     % new snow dendricity 
gspNew = 0.5;   % new snow sphericity 

if P > 0
    % determine initial mass
   mInit = d .* dz;

    % if snow
    if T_air <= 273.15
        z_snow = P/dSnow;               % depth of snow

        % if snow depth is greater than specified min dz, new cell created
        if z_snow > dz_min
            T = [T_air; T];             % new cell T
            dz = [z_snow; dz];          % new cell dz
            d = [dSnow; d];             % new cell d
            W = [0; W];                 % new cell W
            a = [a_SNOW; a];            % new cell a
            re = [reNew; re];           % new cell grain size
            gdn = [gdnNew; gdn];        % new cell grain dendricity
            gsp = [gspNew; gsp];        % new cell grain sphericity

        % if snow depth is less than specified minimum dz snow
        else 
            mass = mInit(1) + P;         % grid cell adjust mass
            dz(1) = dz(1) + P/dSnow;    % adjust grid cell depth      
            d(1) = mass / dz(1);        % adjust grid cell density
            
            % adjust variables as a linearly weighted function of mass
            % adjust temperature (assume P is same temp as air)
            T(1) = (T_air * P + T(1) * mInit(1))/mass;
            
            % adjust a, re, gdn & gsp
            a(1) = (a_SNOW * P + a(1) * mInit(1))/mass;
            re(1) = (reNew * P + re(1) * mInit(1))/mass;
            gdn(1) = (gdnNew * P + gdn(1) * mInit(1))/mass;
            gsp(1) = (gspNew * P + gsp(1) * mInit(1))/mass;
        end
        
    % if rain    
    else
        % rain is added by increasing the mass and temperature of the ice
        % of the top grid cell.  Temperatures are set artifically high to
        % account for the latent heat of fusion.  This is the same as
        % directly adding liquid water to the the snow pack surface but
        % makes the numerics easier.
        
        LF = 0.3345E6;  % latent heat of fusion(J kg-1)
        CI = 2102;      % specific heat capacity of snow/ice (J kg-1 k-1)
        
        % grid cell adjust mass
        mass = mInit(1) + P;

        % adjust temperature
        % liquid: must account for latent heat of fusion
        T(1) = (P *(T_air + LF/CI) + T(1) * mInit(1)) / mass;

        % adjust grid cell density
        d(1) = mass / dz(1);

        % if d > the density of ice, d = dIce
        if d(1) > dIce
            d(1) = dIce;           % adjust d
            dz(1) = mass / d(1);    % dz is adjusted to conserve mass
        end
    end
    
    %% check for conservation of mass
    mass = sum(d .* dz);
    mass_diff = mass - sum(mInit) - P;
    mass_diff = round(mass_diff * 100)/100;
    
    if mass_diff > 0
        error('mass not conserved in accumulation function')
    end
end