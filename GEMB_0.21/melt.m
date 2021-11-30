function [sumM, R, d, T, dz, W, mAdd, a, re, gdn, gsp] = ...
    melt(T, d, dz, W, a, dzMin, zMax, zMin, re, gdn, gsp)
%% MELT ROUTINE

% Description:
% computes the quantity of meltwater due to snow temperature in excess of
% 0 deg C, determines pore water content and adjusts grid spacing
 

%% INITIALIZATION
% specify constants
CtoK = 273.15;  % clecius to Kelvin conversion
CI = 2102;      % specific heat capacity of snow/ice (J kg-1 k-1)
LF = 0.3345E6;  % latent heat of fusion(J kg-1)
dPHC = 830;     % pore hole close off density[kg m-3]
dIce = 910;     % density of ice [kg m-3]

% store initial mass [kg] and energy [J]
m = dz .* d;                    % grid cell mass [kg]
EI = m .* T * CI;               % initial enegy of snow/ice
EW = W .* (LF + CtoK * CI);     % initial enegy of water

mSum0 = sum(W) + sum (m);       % total mass [kg]
sumE0 = sum(EI) + sum(EW);      % total energy [J]

% initialize melt and runoff scalars
R = 0;          % runoff [kg]
sumM = 0;       % total melt [kg]
mAdd = 0;       % mass added/removed to/from base of model [kg]
addE = 0;       % energy added/removed to/from base of model [J]

% calculate temperature excess above 0 °C
exsT = max(0, T - CtoK);        % [K] to [°C]

% new grid point center temperature, T [K]
T = T - exsT;

% specify irreducible water content saturation [fraction]
Swi = 0.07;                     % assumed constant after Colbeck, 1974

%% REFREEZE PORE WATER
% check if any pore water
if sum(W) > 0
%                                                                          disp('PORE WATER REFREEZE')
    % calculate maximum freeze amount, maxF [kg]
    maxF = max(0, -((T - CtoK) .* m * CI) / LF);
    
    % freeze pore water and change snow/ice properties
    dW = min(maxF, W);                              % freeze mass [kg]   
    W = W - dW;                                     % pore water mass [kg]
    m = m + dW;                                     % new mass [kg]
    d = m ./ dz;                                    % density [kg m-3]   
    T = T + (dW.*(LF+(CtoK - T)*CI)./(m.*CI));      % temperature [K]
    
    % if pore water froze in ice then adjust d and dz thickness
    d(d > dIce) = dIce;
    dz = m ./ d;  
end

% squeeze water from snow pack
Wi = (910 - d) .* Swi .* (m ./ d);      % irreducible water content [kg]
exsW = max(0, W - Wi);                  % water "squeezed" from snow [kg]

%% MELT, PERCOLATION AND REFREEZE
 
% run melt algorithm if there is melt water or excess pore water
if (sum(exsT) > 0) || (sum(exsW) > 0)
    %                                                                          disp ('MELT OCCURS')
    % check to see if thermal energy exceeds energy to melt entire cell
    % if so redistribute temperature to lower cells (temperature surplus)
    % (maximum T of snow before entire grid cell melts is a constant
    % LF/CI = 159.1342)
    surpT = max(0, exsT - 159.1342);
    
    if sum(surpT) > 0 % bug fixed 21/07/2016
        
%                                                                          disp('T Surplus')
        % calculate surplus energy
        surpE = surpT .* CI .* m;
        i = 1;
        
        while sum(surpE) > 0
            % use surplus energy to increase the temperature of lower cell
            T(i+1) = surpE(i) / m(i+1) / CI + T(i+1);
            
            exsT(i+1) = max(0, T(i+1) - CtoK) + exsT(i+1);
            T(i+1) = min(CtoK, T(i+1));
            
            surpT(i+1) = max(0, exsT(i+1) - 159.1342);
            surpE(i+1) = surpT(i+1) * CI * m(i+1);
 
            % adjust current cell properties (again 159.1342 is the max T)
            exsT(i) = 159.1342;
            surpE(i) = 0;   
            i = i + 1;
        end
    end

    % convert temperature excess to melt [kg]
    M = exsT .* d .* dz * CI / LF;      % melt
    sumM = sum(M);                      % total melt [kg]
    
    % calculate maximum refreeze amount, maxF [kg]
    maxF = max(0, -((T - CtoK) .* d .* dz * CI)/ LF);
 
    % initialize refreeze, runoff, flxDn and dW vectors [kg]
    F = zeros(length(T),1);
    R = F;
    dW = F;
    flxDn = [0; F];
    
    % determine the deepest gird cell where melt/pore water is generated
    X = find((M > 0 | exsW), 1, 'last');
    X(isempty(X)) = 0;
        
    %% meltwater percolation
    for i = 1:length(T)
        % calculate total melt water entering cell
        inM = M(i)+ flxDn(i);
 
        % break loop if there is no meltwater and if depth is > mw_depth
        if inM == 0 && i > X
            break
 
            % if reaches impermeable ice layer all liquid water runs off (R)
        elseif d(i) >= dIce   % dPHC = pore hole close off [kg m-3]
%                                                                          disp('ICE LAYER')
            % no water freezes in this cell
            % no water percolates to lower cell
            % cell ice temperature & density do not change
            
            m(i) = m(i) - M(i);                     % mass after melt
            Wi = (910-d(i)) * Swi * (m(i)/d(i));    % irreducible water 
            dW(i) = min(inM, Wi - W(i));            % change in pore water
            R(i) = max(0, inM - dW(i));             % runoff
            
            % check if no energy to refreeze meltwater     
        elseif maxF(i) == 0
%                                                                          disp('REFREEZE == 0')
            % no water freezes in this cell
            % cell ice temperature & density do not change
            
            m(i) = m(i) - M(i);                     % mass after melt
            Wi = (910-d(i)) * Swi * (m(i)/d(i));    % irreducible water 
            dW(i) = min(inM, Wi-W(i));              % change in pore water
            flxDn(i+1) = max(0, inM-dW(i));         % meltwater out
            F(i) = 0;                               % no freeze 
 
            % some or all meltwater refreezes
        else
            % change in density density and temperature
%                                                                          disp('MELT REFREEZE')
            %-----------------------melt water-----------------------------
            dz_0 = m(i)/d(i);          
            dMax = (dIce - d(i))*dz_0;              % d max = dIce
            F1 = min([inM; dMax; maxF(i)]);         % maximum refreeze               
            m(i) = m(i) + F1;                       % mass after refreeze
            d(i) = m(i)/dz_0;
            
            %-----------------------pore water-----------------------------
            Wi = (910-d(i))* Swi * dz_0;            % irreducible water 
            dW(i) = min(inM - F1, Wi-W(i));         % change in pore water
            if -dW(i)>W(i)                        
                dW(i) = -W(i);                       % asg changed from dW(i)= W(i) to dW(i)= -W(i) on Nov 14   
            end
            F2 = 0;                                 
            
            %% ---------------- THIS HAS NOT BEEN CHECKED-----------------_
            if dW(i) < 0                            % excess pore water
                dMax = (dIce - d(i))*dz_0;          % maximum refreeze                                             
                maxF2 = min(dMax, maxF(i)-F1);      % maximum refreeze
                F2 = min(-dW(i), maxF2);            % pore water refreeze
                m(i) = m(i) + F2;                   % mass after refreeze
                d(i) = m(i)/dz_0;
            end
            % -------------------------------------------------------------
            
            
            flxDn(i+1) = inM - F1 - dW(i) - F2;     % meltwater out        
            T(i) = T(i) + ...                       % change in temperature
                ((F1+F2)*(LF+(CtoK - T(i))*CI)./(m(i).*CI));
            
            % check if an ice layer forms 
            if d(i) == dIce
%                                                                          disp('ICE LAYER FORMS')
                % excess water runs off
                R(i) = flxDn(i+1);
                % no water percolates to lower cell
                flxDn(i+1) = 0;
            end
        end
    end

%% GRID CELL SPACING AND MODEL DEPTH

if sum(W < 0)
    error ('negative pore water generated in melt equations')
end
   % delete all cells with zero mass
   % adjust pore water
    W = W + dW;   
    
    % delete all cells with zero mass
    D = (m == 0); 
    m(D) = []; W(D) = []; d(D) = []; T(D) = []; a(D) = []; re(D) = []; 
    gdn(D) = []; gsp(D) = [];
 
    % calculate new grid lengths
    dz = m ./ d;
end
 

%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zY = 1.025;
dzMin = zeros(size(dz))+dzMin;
z = cumsum(dz);
n1 = sum(z<10)+1;
n2 = length(z);
dzMin(n1:n2) = zY.^(1:(n2-n1+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if depth is too small
X = find(dz < dzMin,1, 'last');
X(isempty(X)) = 0;

% check to see if last cell it too small 
% Modified to account for mergeing of last cell 08/26/2016
if X == length(dZ)
    lastCellFlag = true;
    X = X-1;
else 
    lastCellFlag = false;
end

%%%
for i = 1:X
    if dz (i) < dzMin(i)                               % merge top cells
        %                                                                          disp('dz > dzMin * 2')
        % adjust variables as a linearly weighted function of mass
        m_new = m(i) + m(i+1);
        T(i+1) = (T(i)*m(i) + T(i+1)*m(i+1)) / m_new;
        a(i+1) = (a(i)*m(i) + a(i+1)*m(i+1)) / m_new;
        re(i+1) = (re(i)*m(i) + re(i+1)*m(i+1)) / m_new;
        gdn(i+1) = (gdn(i)*m(i) + gdn(i+1)*m(i+1)) / m_new;
        gsp(i+1) = (gsp(i)*m(i) + gsp(i+1)*m(i+1)) / m_new;
        
        % merge with underlying grid cell and delete old cell
        dz (i+1) = dz(i) + dz(i+1);                 % combine cell depths
        d(i+1) = m_new / dz(i+1);                   % combine top densities
        W(i+1) = W(i+1) + W(i);                     % combine liquid water
        m(i+1) = m_new;                             % combine top masses
        
        % set cell to 99999 for deletion
        m(i) = 99999;
    end
end

% last cell need to be merged
% Modified to account for mergeing of last cell 08/26/2016

if lastCellFlag
    % find closest cell to merge with
    foo = find(m ~= 99999,2, 'last');
    i = foo(1);
    j = foo(2);
    
    % adjust variables as a linearly weighted function of mass
    m_new = m(i) + m(j);
    T(j) = (T(i)*m(i) + T(j)*m(j)) / m_new;
    a(j) = (a(i)*m(i) + a(j)*m(j)) / m_new;
    re(j) = (re(i)*m(i) + re(j)*m(j)) / m_new;
    gdn(j) = (gdn(i)*m(i) + gdn(j)*m(j)) / m_new;
    gsp(j) = (gsp(i)*m(i) + gsp(j)*m(j)) / m_new;
    
    % merge with underlying grid cell and delete old cell
    dz (j) = dz(i) + dz(j);                 % combine cell depths
    d(j) = m_new / dz(j);                   % combine top densities
    W(j) = W(j) + W(i);                     % combine liquid water
    m(j) = m_new;                             % combine top masses
    
    % set cell to 99999 for deletion
    m(i) = 99999;
end

% delete combined cells
D = (m == 99999);
m(D) = []; W(D) = []; dz(D) = []; d(D) = []; T(D) = []; a(D) = [];
re(D) = []; gdn(D) = []; gsp(D) = [];

% check if any of the top 10 cell depths are too large
X = find(dz(1:10) > 2*dzMin,1, 'last');
X(isempty(X)) = 0;

while i <= X
    if dz(i) > dzMin * 2
%                                                                          disp('dz > dzMin * 2')
       % split in two
        dz = [dz(1:i-1) ; dz(i)/2 ; dz(i)/2 ; dz(i+1:end)];
        W = [W(1:i-1) ; W(i)/2 ; W(i)/2 ; W(i+1:end)]; 
        m = [m(1:i-1) ; m(i)/2 ; m(i)/2 ; m(i+1:end)]; 
        R = [R(1:i-1) ; R(i)/2 ; R(i)/2 ; R(i+1:end)]; 
        T = [T(1:i-1) ; T(i) ; T(i) ; T(i+1:end)];
        d = [d(1:i-1) ; d(i) ; d(i) ; d(i+1:end)]; 
        a = [a(1:i-1) ; a(i) ; a(i) ; a(i+1:end)];
        re = [re(1:i-1) ; re(i) ; re(i) ; re(i+1:end)]; 
        gdn = [gdn(1:i-1) ; gdn(i) ; gdn(i) ; gdn(i+1:end)];
        gsp = [gsp(1:i-1) ; gsp(i) ; gsp(i) ; gsp(i+1:end)]; 
        
        X = X+1;
    else
        i = i+1;
    end
end

 
%% CORRECT FOR TOTAL MODEL DEPTH
% WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
% INTERPRETATION
 
% % calculate total model depth
% z = sum(dz);
% 
% if z < zMin % check if model is too shallow                                                                      
%                                                                          disp('z < zMin')
%     % mass and energy to be added
%     mAdd = m(end) + W(end);
%     addE = T(end) * m(end) * CI;
%     
%     % add a grid cell of the same size and temperature to the bottom
%     dz = [dz; dz(end)];
%     T = [T; T(end)];
%     W = [W; W(end)];
%     m = [m; m(end)];
%     d = [d; d(end)];
%     a = [a; a(end)];
%     re = [re; re(end)];
%     gdn = [gdn; gdn(end)];
%     gsp = [gsp; gsp(end)];
%
% elseif z > zMax % check if model is too deep                                                                                                                                        
%                                                                          disp('z > zMax')
%     % mass and energy loss
%     mAdd = -(m(end) + W(end));
%     addE = -(T(end) * m(end) * CI);
%     
%     % add a grid cell of the same size and temperature to the bottom
%     dz(end) = []; T(end) = []; W(end) = []; m(end) = []; 
%     d(end) = []; a(end) = [];  re(end) = [];  gdn(end) = [];  
%     gsp(end) = [];
% end
 
%% CHECK FOR MASS AND ENERGY CONSERVATION

% calculate final mass [kg] and energy [J]
R = sum(R);
EI = m .* T * CI;
ER = R .* (LF + CtoK * CI);
EW = W .* (LF + CtoK * CI);

mSum1 = sum(W) + sum (m) + R;
sumE1 = sum(EI) + sum(EW);
 
dm = round(mSum0 - mSum1 + mAdd);
dE = round(sumE0 - ...
    sumE1 - ER +  addE);
 
if dm ~= 0 && dE ~= 0
    error('mass and energy are not conserved in melt equations')
elseif dm ~= 0
    error('mass is not conserved in melt equations')
elseif dE ~= 0
    error('energy is not conserved in melt equations')
end

% W = round(W * 10000)/10000;
if sum(W < 0)
    error ('negative pore water generated in melt equations')
end