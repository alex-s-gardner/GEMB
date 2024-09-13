cd GEMB/TEST_GEMBvsISSM
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

cd ../
warning('off','all')
MASTER_RUN
bechmark_out = load('../TEST_DATA/S2A1D2_000001_test.mat');
master_test_out = load('../TEST_DATA/S2A1D2_000001.mat');
if ~isequaln(bechmark_out,master_test_out)
  error('FAILURE in GEMB/MASTER_RUN.m: GEMB output structure differs from ../TEST_DATA/S2A1D2_000001_test.mat . ')
end
