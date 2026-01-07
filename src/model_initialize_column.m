function [T, dz, d, W, re, gdn, gsp, a, a_diffuse] = model_initialize_column(S, LP)

% initialze column variables 
    dz        = grid_initialize(S.column_ztop, S.column_dztop, S.column_zmax, S.column_zy);
    m         = length(dz);
    T         = zeros(m,1) + LP.T_air_mean; % initial grid cell temperature to the annual mean temperature [K]
    d         = zeros(m,1) + S.density_ice; % density to that of ice [kg m-3]
    W         = zeros(m,1);                 % water content to zero [kg m-2]
    re        = zeros(m,1) + 2.5;           % grain size to old snow [mm]
    gdn       = zeros(m,1);                 % grain dentricity to old snow
    gsp       = zeros(m,1);                 % grain sphericity to old snow
    a         = zeros(m,1) + S.albedo_snow; % albedo equal to fresh snow [fraction]
    a_diffuse = zeros(m,1) + S.albedo_snow; % albedo equal to fresh snow [fraction]  
end