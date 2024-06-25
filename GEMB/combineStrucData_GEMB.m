function S = combineStrucData_GEMB(S,LP,runIdx)
% combine input data into structure for passing to parfor loop

% transfer local parameters to S structure
fn = fields(LP);
for i = 1:length(fn)
    S.(fn{i}) = LP.(fn{i});
end

S.runID = [S.runPfx '_' sprintf('%06d', runIdx)];
