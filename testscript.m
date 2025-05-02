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
disp('The thermo function looks good.')

%% Melt test 

load('melt_test_input.mat')
foo_out = load('melt_test_output.mat'); 

[M, Msurf, R, F, T, d, dz, W, mAdd, ~, a, adiff, re, gdn, gsp] = melt(T, d, dz, W, Ra, a, adiff,...
    S.dzMin, S.zMax, S.zMin, S.zTop, S.zY, re, gdn, gsp, dIce);

assert(isequal(M,foo_out.M),'M might be wrong.')
assert(isequal(Msurf,foo_out.Msurf),'M might be wrong.')
assert(isequal(R,foo_out.R),'R might be wrong.')
assert(isequal(F,foo_out.F),'F might be wrong.')
assert(isequal(T,foo_out.T),'T might be wrong.')
assert(isequal(d,foo_out.d),'d might be wrong.')
assert(isequal(dz,foo_out.dz),'dz might be wrong.')
assert(isequal(W,foo_out.W),'W might be wrong.')
assert(isequal(mAdd,foo_out.mAdd),'mAdd might be wrong.')
assert(isequal(a,foo_out.a),'a might be wrong.')
assert(isequal(adiff,foo_out.adiff),'adiff might be wrong.')
assert(isequal(dz,foo_out.dz),'dz might be wrong.')
assert(isequal(re,foo_out.re),'re might be wrong.')
assert(isequal(gdn,foo_out.gdn),'gdn might be wrong.')
assert(isequal(gsp,foo_out.gsp),'gsp might be wrong.')
disp('The melt function looks good.')

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
