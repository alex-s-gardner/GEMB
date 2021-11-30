function dz = gridInitialize(zTop, dzTop, sMax, zY)

% Description: 
% This file sets up the initial grid spacing and total grid depth.  The
% grid structure is set as constant grid length 'dzTop' for the top
% 'zTop' meters of the model grid. Bellow 'zTop' the gid length increases
% linearly with depth

%----------------------Calculate Grid Lengths------------------------------
% calculate number of top grid points
gpTop = zTop/dzTop;

% check to see if the top grid cell structure length (dzTop) goes evenly 
% into specified top structure depth (zTop)
if gpTop ~= round(gpTop)
    error(['top grid cell structure length does not go evenly into ' ...
        'specified top structure depth, adjust dzTop or zTop'])

% make sure top grid cell structure length (dzTop) is greater than 5 cm
elseif dzTop < 0.05
    warning('initial top grid cell length (dzTop) is < 0.05 m')
end

% initialize top grid depth vector
dzT = ones(gpTop,1)*dzTop;

% bottom grid cell depth = x*zY^(cells from to structure)

% initialize bottom vectors
dzB = zeros(((sMax - zTop)/dzTop),1);
gp0 = dzTop;
z0 = zTop;
i = 1;

while sMax > z0
    dzB(i,1) = gp0 * zY;
    gp0 = dzB(i,1);
    z0 = z0 + gp0;
    i = i + 1;
end

% delete excess cells from bottom vector 
dzB(dzB == 0) = [];

% combine top and bottom dz vectors
dz = [dzT ; dzB];

%% ---------NEED TO IMPLEMENT A PROPER GRID STRECHING ALGORITHM------------
