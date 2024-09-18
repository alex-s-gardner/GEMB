% This script loads some input test data, feeds it to thermo.m, then 
% verifies that the results are what's expected.
% Chad Greene, NASA/JPL, September 2024. 

%% Load test data

load('thermo_test_input.mat')

foo = load('thermo_test_output.mat'); 

%% Call thermo.m 

[shf_cum, lhf_cum, T, EC, ulwrf] = thermo(T, re, dz, d, swf, dlwrf, Ta, V, eAir, pAir, tcIdx, eIdx, ...
        teValue, dulwrfValue, teThresh, Ws, dt0, dzMin, Vz, Tz, dtScaling, dIce, isdeltaLWup);

%% Check all the outputs: 

assert(isequal(shf_cum,foo.shf_cum),'shf_cum might be wrong.')
assert(isequal(lhf_cum,foo.lhf_cum),'lhf_cum might be wrong.')
assert(isequal(T,foo.T),'T might be wrong.')
assert(isequal(EC,foo.EC),'EC might be wrong.')
assert(isequal(ulwrf,foo.ulwrf),'ulwrf might be wrong.')
disp('If the code made it this far, everything looks good!')