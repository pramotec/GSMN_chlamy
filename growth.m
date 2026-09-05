%% ========== Supplementary Growth Analysis ==========
% Original GSMN + biomass objective + light upper bound

model_orig = model_final;

% Check biomass reaction
biomassRxn = 'Biomass_Chlamy_mixo';

if ~ismember(biomassRxn, model_orig.rxns)
    error('%s is not present in the original model.', biomassRxn);
end

% Check light reaction
lightRxn = 'EX_photonVis_e';

if ~ismember(lightRxn, model_orig.rxns)
    error('%s is not present in the original model.', lightRxn);
end

% Clean infinite bounds
model_orig.lb(isinf(model_orig.lb) & model_orig.lb < 0) = -1000;
model_orig.ub(isinf(model_orig.ub) & model_orig.ub > 0) = 1000;

% Set biomass as objective
model_orig = changeObjective(model_orig, biomassRxn);

%% HL
model_HL_orig = changeRxnBounds(model_orig, lightRxn, -1000, 'l');
model_HL_orig = changeRxnBounds(model_HL_orig, lightRxn, 0, 'u');

%% LL
model_LL_orig = changeRxnBounds(model_orig, lightRxn, -50, 'l');
model_LL_orig = changeRxnBounds(model_LL_orig, lightRxn, 0, 'u');

%% FBA
sol_HL_orig = optimizeCbModel(model_HL_orig, 'max');
sol_LL_orig = optimizeCbModel(model_LL_orig, 'max');

fprintf('\n========== Supplementary Growth Analysis ==========\n');
fprintf('Objective: %s\n', biomassRxn);
fprintf('HL photon upper capacity: 1000\n');
fprintf('LL photon upper capacity: 50\n\n');

fprintf('HL theoretical maximum growth rate: %.6f\n', sol_HL_orig.f);
fprintf('LL theoretical maximum growth rate: %.6f\n', sol_LL_orig.f);

if sol_LL_orig.f ~= 0
    fprintf('HL/LL growth-rate ratio: %.3f\n', ...
        sol_HL_orig.f / sol_LL_orig.f);
end

fprintf('===================================================\n');