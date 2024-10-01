cd GEMB
warning('off','all')

%% thermo test 

load('../TEST_DATA/thermo_test_input.mat')

foo_out = load('../TEST_DATA/thermo_test_output.mat'); 

% Call thermo.m 
[shf_cum, lhf_cum, T, EC, ulwrf] = thermo(T, re, dz, d, swf, dlwrf, Ta, V, eAir, pAir, tcIdx, eIdx, ...
        teValue, dulwrfValue, teThresh, Ws, dt0, dzMin, Vz, Tz, dtScaling, dIce, isdeltaLWup);

% Check all the outputs: 
assert(isequal(shf_cum,foo_out.shf_cum),'shf_cum might be wrong.')
assert(isequal(lhf_cum,foo_out.lhf_cum),'lhf_cum might be wrong.')
assert(isequal(T,foo_out.T),'T might be wrong.')
assert(isequal(EC,foo_out.EC),'EC might be wrong.')
assert(isequal(ulwrf,foo_out.ulwrf),'ulwrf might be wrong.')

%%
MASTER_RUN
bechmark_out = load('../TEST_DATA/S2A1D2_000001_test.mat');
master_test_out = load('../TEST_DATA/S2A1D2_000001.mat');

cd TEST_GEMBvsISSM
testGEMBISSM
if max_meltbias_with_gemb_matlab
  error('FAILURE in GEMB/TEST_GEMBvsISSM/testGEMBISSM.m: new melt values do not match archive GEMBtest_output.mat . ')
end
if final_layerbias_with_gemb_matlab
  error('FAILURE in GEMB/TEST_GEMBvsISSM/testGEMBISSM.m: the run ends with numbers of layers that do not match archive GEMBtest_output.mat. ')
end
if ~isequaln(a,b)
  error('FAILURE in GEMB/TEST_GEMBvsISSM/testGEMBISSM.m: GEMB output structure differs from GEMBtest_output.mat . ')
end
if ~isequaln(bechmark_out,master_test_out)
  error('FAILURE in GEMB/MASTER_RUN.m: GEMB output structure differs from ../TEST_DATA/S2A1D2_000001_test.mat . ')
end
