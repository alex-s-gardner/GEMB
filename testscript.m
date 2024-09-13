cd GEMB/TEST_GEMBvsISSM
testGEMBISSM
if max_meltbias_with_gemb_matlab
  error('FAILURE in GEMB/TEST_GEMBvsISSM/testGEMBISSM.m: new melt values do not match archive GEMBtest_output.mat . ')
end
if final_layerbias_with_gemb_matlab
  error('FAILURE in GEMB/TEST_GEMBvsISSM/testGEMBISSM.m: the run ends with numbers of layers that do not match ...
    archive GEMBtest_output.mat. ')
end
