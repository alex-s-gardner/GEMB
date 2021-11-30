% number of processors to use
mod.NW = 48;
mod.PC = 'off'; %'flux'; %'flux'; % 'on';
TEST = true;

%% GEMB INITIALIZATION
% unique model run ID to save output as
S.runPfx = 'S2A1D2';

% spin-up
S.spinUp = 2;   % number of cycles of met data run before output
% calcualted, set spinUp = 0 for no spin up

% select method of calculating albedo and subsurface absorption
%   1 : effective grain radius (Gardner & Sharp, 2009)
%   2 : effective grain radius (Brun et al., 2009)
%   3 : density and cloud amount (Greuell & Konzelmann, 1994)
%   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)
S.aIdx = 1;

% apply all SW to top grid cell (0) or allow SW to penetrate surface (1)
S.swIdx = 1;

% select densification model to use:
%       1 = emperical model of Herron and Langway (1980)
%       2 = semi-emerical model of Anthern et al. (2010)
%       3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
%       4 = DO NOT USE: emperical model of Li and Zwally (2004)
%       5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
S.denIdx = 2;

% GRID INITIALIZATION
% set depth of top grid cell structure (constant grid length) [m]
S.zTop = 10;

% set initial top vertical and min allowable grid spacings [m]
S.dzTop = 0.05; %0.05;
S.dzMin = S.dzTop / 2;

% set initial/max and min model depth [m]
S.zMax = 250;
S.zMin = ceil(S.zMax/2 /10)*10;

% strech grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]
S.zY = 1.10;

% OTHER
% specify frequency to output data (density, grid length, and temperature)
%   - 'monthly'
%   - 'daily'
S.outputFreq = 'monthly';

% input data directory
if TEST
    S.inputDIR = fullfile('..','TEST/');
    S.outDIR = fullfile('..','TEST/');
else
    S.inputDIR = '/Volumes/MasterBrain/data/GEMB/CFSR/T62/';
    % output directory
    S.outDIR = '/Volumes/MasterBrain/data/GEMB/Output/';
end

% FIXED ALBEDO VARIABLES
switch S.aIdx
    case {1 2}
        % albedo tuning parameters
        % for methods of calculating albedo see albedo function
        % METHOD 1 & 2
        S.aSnow = 0.85;         % new snow albedo (0.64 - 0.89)
        S.aIce = 0.48;          % range 0.27-0.58 for old snow
        
    case 3
        % RADIATION CORRECTION FACTORS
        % -> only used for met station data and Greuell & Konzelmann, 1994 albedo
        S.cldFrac = 0.1;        % average cloud amount
        
    case 4
        % additonal tuning parameters albedo as a funtion of age and water content
        % (Bougamont et al., 2005)
        S.t0wet = 15;           % time scale for wet snow (15-21.9) [d]
        S.t0dry = 30;           % warm snow timescale (30) [d]
        S.K = 7;                % time scale temperature coef. (7) [d]
end

%% RUN GEMB
runPfx = S.runPfx;
inputDIR = S.inputDIR;

% open matlab pool for parallel processing
if strcmp(mod.PC, 'on')
    parpool(mod.NW)
elseif strcmp(mod.PC, 'flux')
    parpool('mpiexecflux', mod.NW)
end

if TEST
    I = load(fullfile(S.inputDIR,'TEST_INPUT_1'));
    S0 = combineStrucData_GEMB(S,I.LP,1);
    GEMB(I.P0,I.Ta0,I.V0,I.dateN,I.dlw0,I.dsw0,I.eAir0,I.pAir0, S0)
else
    H = dir(fullfile(S.inputDIR,'input*'));
    parfor  runIdx = 1:length(H)
        I = load(fullfile(S.inputDIR,['input_' sprintf('%06d', runIdx)]));
        S0 = combineStrucData_GEMB(S,I.LP,runIdx);
        GEMB(I.P0,I.Ta0,I.V0,I.dateN,I.dlw0,I.dsw0,I.eAir0,I.pAir0, S0)
    end
    
    delete(gcp)
end