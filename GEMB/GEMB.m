function GEMB(P0, Ta0, V0, dateN, dlw0, dsw0, eAir0, pAir0, S, isrestart)
%% Glacier Energy and Mass Balance (GEMB) Model

%% Model Documentation
% *Created by*: Alex S. Gardner, JPL/Caltech
% 
% Date last modified: May 19, 2015                                    
%%
% *Description:*                                                             
%%
% This program calculates a 1-D surface glacier mass balance and 
% includes detailed representation of subsurface process. Key features:
%
% * melt water percolation and refreeze
% * pore water retention                                               
% * dynamic albedo with long-term memory                               
% * subsurface temperature diffusion                                   
% * subsurface penetration of shortwave radiation                       
%%
%
% *Reference:*
%
% Many aspects of the general model structure and specific 
% parameterizations are taken from:
%
% Bassford, R. P., 2002: Geophysical and numerical modelling investigations
% of the ice caps on Severnaya Zemlya, Bristol Glaciology Centre, School of
% Geographical Sciences, University of Bristol 220.
%
% Bougamont, M. and J. L. Bamber, 2005: A surface mass balance model for
% the Greenland Ice Sheet. Journal of Geophysical Research-Earth Surface,
% 110, doi:10.1029/2003JD004451.
%
% Greuell, W.; Konzelmann, T., 1994. Numerical Modeling of the
% Energy-Balance and the Englacial Temperature of the Greenland Ice-Sheet -
% Calculations for the Eth-Camp Location (west Greenland, 1155M Asl).
% Global and Planetary Change. 9(1994)
%%
%
% *Model variables:*
%
% * z   = grid cell depth below surface [m]
% * dz  = grid cell size [m]
% * a   = albedo [fraction]
% * EC  = evaporation (-) & condensation (+) [kg m^-^2]
% * d   = denity [kg m^-^3]
% * mass   = mass [kg m^-^2]
% * M   = melt [kg m^-^2]
% * W   = liquid water content [kg m^-^2]
% * re  = optically equivelant grain radius [mm]
% * gdn = grain dentricity
% * gsp = grain sphericity
% * R   = runoff [kg m-2] 
% * ps  = pore space [m^3/m^2] 
%
% *Output variables:*
%
% * comp1 = dry snow ccompaction [m]
% * comp2 = melt compaction [m]
%
%% Input
%
% * dateN : date/time [UTC: days since ...]
% * Ta: 2m air temperature [K]
% * V: wind speed [m s^-^1]
% * dsw: downward shortwave radiation flux [W m^-^2]
% * dlw: downward longwave radiation flux [W m^-^2] 
% * P: precipitation [kg m^-^2]
% * eAir: screen level vapor pressure [Pa]
% * pAir: surface pressure [Pa]
% * S.Tmean: mean annual temperature [K]
% * S.C = mean annual snow accumulation [kg m-2 yr-1]
% * S.Tz: height above ground at which temperature (T) was sampled [m]
% * S.Vz: height above ground at which wind (V) eas sampled [m]

% specify met station 
    % * automatically generated met data = 0
    % * Devon Site 1,2,3 = 1,2,3
    % * CFSR extracted climate forcing = 999
% S.site = 999; 

% Daily (1) or hourly (2) data
% S.dataFreq = 2;
% if S.site < 999
%      S.saveFreq = 24
%     [dateN, Ta0, V0, dsw0, dlw0, P0, eAir0, pAir0, dt, LP] = ...
%         loadMetData(S.cldFrac, S.site, S.dataFreq, S.saveFreq);
% else
%     runIdx = 1;
%     load(fullfile(S.cfsrForcingDIR,['input_' sprintf('%06d', runIdx)]))
% end

disp(['------------------ STARTING RUN # ' num2str(S.runID) ' --------------------' ])
tic                                     % start timer
dt = (dateN(2)-dateN(1)) * (60*60*24);  % input time step in seconds

% test switches
checkInput      = false;
overridePrecip  = false;

% Constants
dIce = 910;     % density of ice [kg m-3]

if overridePrecip
    % override pricipitiation data for testing
    disp('--> setting precipitation to annual mean <--')
    P0(:) = S.C*(dt/(60*60*24*364.25)); % FIX NUM DAYS IN YEAR
end

if checkInput
    % display input parameter ranges to screen for sanity check
    disp(' - - - - - check input forcing - - - - - ')
    X = {'dateN', 'Ta0', 'V0', 'dsw0', 'dlw0', 'P0', 'eAir0', 'pAir0', 'dt'};
    for i = 1:length(X)
        
        fprintf('%s \t : %6.2f - %6.2f \t size : [%7.0f %7.0f]\n', ...
            X{i}, min(eval(X{i})), max(eval(X{i})), size(eval(X{i}),1), ...
            size(eval(X{i}),2))
    end
end

%% Generate model grid

dz = gridInitialize(S.zTop, S.dzTop, S.zMax, S.zY);

%% Initialize model variables

% --------------- INPORTANT NOTE ABOUT GRAIN PROPERTIES -------------------
% initial grain properties must be chosen carefully since snow with a 
% density that exceeds 400 kg m-3 will no longer undergo metamorphosis and
% therefore if grain properties are set inappropriately they will be
% carried all through the model run.  
%
% !!!! grainGrowth model needs to be fixed to allow evolution of snow !!!!!
% !!!!  grains for densities > 400 kg m-3                             !!!!!  
% -------------------------------------------------------------------------

% initialize profile variables
if ~isrestart
    m = length(dz);
    Z = zeros(m,1);         % create zeros matrix
    d = Z + dIce;           % density to that of ice [kg m-3]
    re = Z + 2.5;           % set grain size to old snow [mm]
    gdn = Z;                % set grain dentricity to old snow 
    gsp = Z;                % set grain sphericity to old snow 
    EC = 0;                 % surface evaporation (-) condensation (+) [kg m-2]
    W = Z;                  % set water content to zero [kg m-2]
    M = Z;                  % set melt water to zero [kg m-2]
    Ra = Z;                 % set rain amount to zero [kg m-2]
    F = Z;                  % set refreeze to zero [kg m-2]
    a = Z + S.aSnow;        % set albedo equal to fresh snow [fraction]
    adiff = Z + S.aSnow;    % set albedo equal to fresh snow [fraction]
    Msurf = 0.0;            % initialize surface melt for albedo parameterization

    % set initial grid cell temperature to the annual mean temperature [K]
    T = Z + S.Tmean;
else
    m = S.Sizeini;
    Z = zeros(m,1);         % create zeros matrix
    dz = S.Dzini;           % set layering
    d = S.Dini;             % set density [kg m-3]
    re = S.Reini;           % set grain size [mm]
    gdn = S.Gdnini;         % set grain dentricity
    gsp = S.Gspini;         % set grain sphericity
    EC = S.ECini;           % surface evaporation (-) condensation (+) [kg m-2]
    W = S.Wini;             % set water content [kg m-2]
    T = S.Tini;             % set snow temperature [K]
    M = Z;                  % set melt water to zero [kg m-2]
    Ra = Z;                 % set rain amount to zero [kg m-2]
    F = Z;                  % set refreeze to zero [kg m-2]
    a = S.Aini;             % set albedo [fraction]
    adiff = S.Adiffini;     % set albedo [fraction]
    Msurf = 0.0;            % initialize surface melt for albedo parameterization
end

% fixed lower temperature bounday condition - T is fixed
T_bottom = T(end);  

% deteremine save time steps
dateV = datevec([dateN; (dateN(end) + dateN(end)-dateN(end-1))]);
switch S.outputFreq
    case 'monthly'
        outIdx = (dateV(1:end-1,2) - dateV(2:end,2)) ~= 0;
    case 'daily'
        outIdx = (dateV(1:end-1,3) - dateV(2:end,3)) ~= 0;
    case '3hourly'
        outIdx = (dateV(1:end-1,4) - dateV(2:end,4)) ~= 0;
end

% initialize output structure

% single level time series
S.varName.monolevel = {'time', 'Ta', 'P', 'M', 'R', 'F', 'EC', 'netSW', ...
    'netLW', 'shf', 'lhf', 'a1', 'netQ', 're1', 'd1', 'm', 'FAC'}; 

Z       = nan(1,sum(outIdx));    
O.time  = dateN(outIdx)';
O.M     = Z; 
O.R     = Z; 
O.F     = Z; 
O.netSW = Z; 
O.netLW = Z; 
O.shf   = Z; 
O.lhf   = Z; 
O.a1    = Z; 
O.netQ  = Z; 
O.re1   = Z; 
O.d1    = Z;
O.Ta    = Z; 
O.P     = Z; 
O.comp1 = Z; 
O.comp2 = Z; 
O.ps    = Z;
O.m     = Z;

% time averages/totals
I = find(outIdx);                      % save index
for i = 1:length(I)
    if i == 1
        O.Ta(i) = mean(Ta0(1:I(i))) - 273.15;  % convert from K to deg C
        O.P(i)  = sum(P0(1:I(i)));
    else
        O.Ta(i) = mean(Ta0((I(i-1)+1):I(i))) - 273.15;
        O.P(i)  = sum(P0((I(i-1)+1):I(i)));
    end
end

% multi level time series
S.varName.profile = {'d', 'T', 'W', 'a', 'dz', 're', 'gdn', 'gsp'};
S.addCells        = 10000;   % number of addtional vertical levels

Z     = nan(length(d)+S.addCells,length(I));       
O.d   = Z; 
O.T   = Z; 
O.W   = Z; 
O.dz  = Z; 
O.re  = Z;
O.gdn = Z; 
O.gsp = Z; 
O.ps  = Z;

% clear unwanted variables
clear X Z n dateV fn

% initialize cumulative output values
OV.varName = {'R', 'M', 'F', 'P', 'EC', 'Ra', 'mAdd', 'netSW', 'netLW', 'shf', ...
    'lhf', 'a1', 're1', 'ulw', 'd1', 'comp1', 'comp2', 'm', 'netQ', 'FAC'};

% set cumulative values zero
for v = 1:length(OV.varName)
    OV.(OV.varName{v}) = 0;
end

OV.count = 0;

%% Start year loop for model spin up
for yIdx = 1:S.spinUp + 1

    % Determine initial mass [kg]:
    initMass = sum (dz .* d) + sum(W);
    
    % Initialize cumulative variables:
    sumR       = 0; 
    sumF       = 0; 
    sumM       = 0; 
    sumEC      = 0; 
    sumP       = 0; 
    sumMassAdd = 0; 
    sumMsurf   = 0; 
    sumRa      = 0;

    %% Start loop for data frequency
    
    % Specify the time range over which the mass balance is to be calculated:
    for dIdx = 1:length(dateN)
        
        % Extract daily data:
        dlw  = dlw0(dIdx);     % downward longwave radiation flux [W m-2]
        dsw  = dsw0(dIdx);     % downward shortwave radiation flux [W m-2]
        Ta   = Ta0(dIdx);      % screen level air temperature [K]
        P    = P0(dIdx);       % precipitation [kg m-2]
        V    = V0(dIdx);       % wind speed [m s-1]
        eAir = eAir0(dIdx);    % screen level vapor pressure [Pa]
        pAir = pAir0(dIdx);    % screen level air pressure [Pa]

        % albedo calculations contained in switch to minimize passing of
        % variables to albedo function
        switch S.aIdx
            case {1,2}
                % if we are provided with cc and cot values, extract for the timestep
                if numel(S.ccsnowValue)>1
                    ccsnowValue = S.ccsnowValue(dIdx);
                else
                    ccsnowValue = S.ccsnowValue;
                end

                if numel(S.cciceValue)>1
                    cciceValue = S.cciceValue(dIdx);
                else
                    cciceValue = S.cciceValue;
                end

                if numel(S.cotValue)>1
                    cotValue = S.cotValue(dIdx);
                else
                    cotValue = S.cotValue;
                end

                if numel(S.szaValue)>1
                    szaValue = S.szaValue(dIdx);
                else
                    szaValue = S.szaValue;
                end

                if numel(S.dswdiffrf)>1
                    dswdiffrf = S.dswdiffrf(dIdx);
                else
                    dswdiffrf = S.dswdiffrf;
                end

                % snow grain metamorphism
                [re, gdn, gsp]  = ...
                    grainGrowth(T, dz, d, W, re, gdn, gsp, dt, S.aIdx);

                % calculate snow, firn and ice albedo
                [a, adiff] = albedo(S.aIdx, re, dz, d, [], S.aIce, S.aSnow, S.aValue, S.adThresh, a, adiff, T, W, P, EC, ...
                    Msurf, ccsnowValue, cciceValue, szaValue, cotValue, [], [], [], dt, dIce);

            case 3   
                if numel(S.cldFrac)>1
                    cldFrac = S.cldFrac(dIdx);
                else
                    cldFrac = S.cldFrac;
                end

                % calculate snow, firn and ice albedo
                [a, adiff] = albedo(S.aIdx, re, dz, d, cldFrac, S.aIce, S.aSnow, S.aValue, S.adThresh,...
                    a, adiff, [], [], [], [], [], [], [], [], [], [], [], [], [], dIce);
                
            case 4
                % calculate snow, firn and ice albedo
                [a, adiff] = albedo(S.aIdx, [], [], d, [], S.aIce, S.aSnow, S.aValue, S.adThresh, a, adiff, T, ...
                    W, P, EC, [], [], [], [], [], S.t0wet, S.t0dry, S.K, dt, dIce);
        end
        
        % determine distribution of absorbed sw radation with depth
        swf = shortwave(S.swIdx, S.aIdx, dsw, dswdiffrf, a(1), adiff(1), d, dz, re, dIce);

        % calculate net shortwave [W m-2]
        netSW = sum(swf);

        % calculate new temperature-depth  profile   
        % and calculate turbulent heat fluxes [W m-2]
        [shf, lhf, T, EC, ulw] = thermo(T, re, dz, d, swf, dlw, Ta, V, eAir, pAir, S.tcIdx, S.eIdx, ...
            S.teValue, S.dulwrfValue, S.teThresh, W(1), dt, S.dzMin, S.Vz, S.Tz, S.ThermoDeltaTScaling, dIce, ...
            S.isdeltaLWup);     

        % change in thickness of top cell due to evaporation/condensation
        % assuming same density as top cell
        % ## NEED TO FIX THIS IN CASE ALL OR MORE OF CELL EVAPORATES ##
        dz(1) = dz(1) + EC / d(1);

        % add snow/rain to top grid cell adjusting cell depth, temperature 
        % and density     
        [T, dz, d, Ra, W, a, adiff, re, gdn, gsp] = accumulation(S.aIdx, S.dsnowIdx, S.Tmean, Ta, T, dz, d, ...
            P, W, S.dzMin, S.C, V, S.Vmean, a, adiff, S.aSnow, re, gdn, gsp, dIce);

        % calculate water production, M [kg m-2] resulting from snow/ice
        % temperature exceeding 273.15 deg K (> 0 deg C), runoff R [kg m-2] 
        % and resulting changes in density and determine wet compaction [m]
        comp2 = sum(dz); 
        [M, Msurf, R, F, T, d, dz, W, mAdd, ~, a, adiff, re, gdn, gsp] = melt(T, d, dz, W, Ra, a, adiff,...
            S.dzMin, S.zMax, S.zMin, S.zTop, S.zY, re, gdn, gsp, dIce);

        comp2 = (comp2 - sum(dz));
        
        % allow non-melt densification and determine compaction [m]
        comp1 = sum(dz); 
        [d, dz] = densification(S.denIdx, S.aIdx, S.swIdx, S.adThresh, d, T, dz, S.C, dt, re, S.Tmean, dIce);
        comp1 = (comp1 - sum(dz));

        % calculate upward longwave radiation flux [W m-2]
        % not used in energy balance 
        % CALCULATED FOR EVERY SUB-TIME STEP IN THERMO EQUATIONS
        % ulw = 5.67E-8 * T(1)^4;

        % calculate net longwave [W m-2]
        netLW = dlw - ulw;

        % sum component mass changes [kg m-2]
        sumMassAdd = mAdd + sumMassAdd;
        sumM       = M + sumM;
        sumMsurf   = Msurf + sumMsurf;
        sumR       = R + sumR;
        sumW       = sum(W);
        sumP       = P +  sumP;
        sumEC      = EC + sumEC;   % evap (-)/cond(+)
        sumRa      = Ra + sumRa;
        sumF       = F + sumF;

        % calculate total system mass
        sumMass = sum(dz .* d);
        FAC     = sum(dz.*(dIce - min(d,dIce)))/1000;
        dMass   = sumMass + sumR + sumW - sumP - sumEC - initMass - sumMassAdd;
        dMass   = round(dMass * 100)/100;
        
        % check mass conservation
        if dMass ~= 0
            error('total system mass not conserved in MB function')
        end

        % check bottom grid cell T is unchanged
        if (T(end)-T_bottom)>1e-8
            warning('T(end)~=T_bottom')
        end

        if yIdx == S.spinUp + 1
            % initialize cumulative and average variables for output
            d1   = d(1); 
            a1   = a(1); 
            re1  = re(1); 
            netQ = netSW + netLW + shf + lhf;

            for v = 1:length(OV.varName)
                OV.(OV.varName{v}) = eval(OV.varName{v}) + OV.(OV.varName{v});
            end

            OV.count = OV.count + 1;

            if outIdx(dIdx)
                % Store model output in structure format

                % time averaged monolevel values
                r = sum(outIdx(1:dIdx));

                for v = 1:length(OV.varName)
                    % check if output is a cumulative value
                    if sum(strcmp(OV.varName{v}, {'M', 'R', 'F', 'EC', 'P', 'Ra','mAdd', 'comp1', 'comp2'})) == 1
                        O.(OV.varName{v})(r) = OV.(OV.varName{v});
                    else
                        % if not cumulative then divide by time steps
                        O.(OV.varName{v})(r) = OV.(OV.varName{v}) / OV.count;
                    end
                end

                % instantaneous level data
                o = (size(d,1) - 1);
                O.re(end-o:end,r)  = re;     
                O.d(end-o:end,r)   = d;
                O.T(end-o:end,r)   = T;       
                O.W(end-o:end,r)   = W;
                O.dz(end-o:end,r)  = dz;     
                O.gdn(end-o:end,r) = gdn;
                O.gsp(end-o:end,r) = gsp;   
                O.ps(end-o:end,r)  = sum(dz) - sumMass/910;
                O.m(r)             = o+1;

                % set cumulative values back to zero
                for v = 1:length(OV.varName)
                    OV.(OV.varName{v}) = 0;
                end
                
                OV.count = 0;
            end
        end
    end

    % display cycle completed and time to screen
    disp([num2str(S.runID) ': cycle: ' num2str(yIdx) ' of '  ...
        num2str(S.spinUp + 1) ', cpu time: ' num2str(round(toc)) ' sec,'...
        ' avg melt: ' num2str(round(sumM/(dateN(end)-dateN(1))*365.25)) ...
        ' kg/m2/yr']);
end

%% Save model output and model run settings
save(fullfile(S.outDIR, S.runID), '-struct', 'O', '-v7.3')
save(fullfile(S.outDIR, S.runID), 'S', '-append')

end