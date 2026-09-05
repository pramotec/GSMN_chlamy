%% SUPPLEMENTARY ANALYSIS
% This analysis uses the original GSMN (model_final)

clearvars -except model_final

%% Identify potential limiting exchange reactions
% Original GSMN without proteomics constraints

testModel = model_final;

testModel = changeRxnBounds( ...
    testModel, 'EX_photonVis_e', -1000, 'l');

testModel = changeObjective( ...
    testModel, 'Biomass_Chlamy_mixo');

sol_base = optimizeCbModel(testModel);

baselineGrowth = sol_base.f;

fprintf('\nBaseline maximum biomass = %.8f\n\n', ...
    baselineGrowth);

exchangeRxns = findExcRxns(testModel);
exRxns = testModel.rxns(exchangeRxns);

results = cell(length(exRxns),5);

for i = 1:length(exRxns)

    rxnID = exRxns{i};
    idx = find(strcmp(testModel.rxns, rxnID));

    originalLB = testModel.lb(idx);

    % Only test uptake reactions
    if originalLB < -1e-9

        perturbedModel = testModel;

        % Double uptake capacity
        perturbedModel.lb(idx) = 2 * originalLB;

        sol = optimizeCbModel(perturbedModel);

        if sol.stat == 1
            newGrowth = sol.f;
            growthIncrease = newGrowth - baselineGrowth;
        else
            newGrowth = NaN;
            growthIncrease = NaN;
        end

    else

        newGrowth = baselineGrowth;
        growthIncrease = 0;

    end

    results{i,1} = rxnID;
    results{i,2} = originalLB;
    results{i,3} = 2 * originalLB;
    results{i,4} = newGrowth;
    results{i,5} = growthIncrease;

end

T_limiting = cell2table(results, ...
    'VariableNames', ...
    {'ExchangeRxn', ...
     'OriginalLB', ...
     'RelaxedLB', ...
     'GrowthAfterRelaxation', ...
     'GrowthIncrease'});

% Sort by growth increase
T_limiting = sortrows( ...
    T_limiting, 'GrowthIncrease', 'descend');

% Show top 10
disp(T_limiting(1:min(10,height(T_limiting)),:));

%% 1. Define model and reactions

model_growth = model_final;

biomassRxn = 'Biomass_Chlamy_mixo';
photonRxn  = 'EX_photonVis_e';
acetateRxn = 'EX_ac_e';
co2Rxn     = 'EX_co2_e';


%% 2. Check reaction IDs

photonIdx  = findRxnIDs(model_growth, photonRxn);
acetateIdx = findRxnIDs(model_growth, acetateRxn);
co2Idx     = findRxnIDs(model_growth, co2Rxn);
biomassIdx = findRxnIDs(model_growth, biomassRxn);

if photonIdx == 0
    error('Photon reaction not found: %s', photonRxn);
end

if acetateIdx == 0
    error('Acetate reaction not found: %s', acetateRxn);
end

if co2Idx == 0
    error('CO2 reaction not found: %s', co2Rxn);
end

if biomassIdx == 0
    error('Biomass reaction not found: %s', biomassRxn);
end


%% 3. Set photon availability

% Photon sensitivity analysis showed that biomass reaches a plateau
% at approximately -67.
%
% Therefore, -80 is used here to provide sufficient photon availability
% without making photon the major limiting constraint.

model_growth.lb(photonIdx) = -80;


%% 4. Set biomass reaction as objective

model_growth = changeObjective(model_growth, biomassRxn);


%% 5. Baseline biomass production

sol_base = optimizeCbModel(model_growth);

if sol_base.stat ~= 1
    error('Baseline biomass optimization failed.');
end

baselineGrowth = sol_base.f;

fprintf('\n============================================\n');
fprintf('Baseline biomass production\n');
fprintf('============================================\n');

fprintf('Photon reaction : %s\n', photonRxn);
fprintf('Photon LB       : %.2f\n', model_growth.lb(photonIdx));

fprintf('Acetate reaction: %s\n', acetateRxn);
fprintf('Acetate LB      : %.2f\n', model_growth.lb(acetateIdx));

fprintf('CO2 reaction    : %s\n', co2Rxn);
fprintf('CO2 LB          : %.2f\n', model_growth.lb(co2Idx));

fprintf('Biomass reaction: %s\n', biomassRxn);
fprintf('Baseline biomass: %.6f\n', baselineGrowth);


%% 6. Define acetate and CO2 uptake capacities

% Negative values represent uptake for exchange reactions.

acetateLevels = [-2 -2.5 -3 -3.5 -4 -4.5];

co2Levels = [-2 -2.5 -3 -3.5 -4 -4.5];


%% 7. Two-dimensional sensitivity analysis

growthMatrix = nan(length(acetateLevels), ...
                   length(co2Levels));

for i = 1:length(acetateLevels)

    for j = 1:length(co2Levels)

        testModel = model_growth;

        % Set acetate uptake
        testModel.lb(acetateIdx) = acetateLevels(i);

        % Set CO2 uptake
        testModel.lb(co2Idx) = co2Levels(j);

        % Optimize biomass
        sol = optimizeCbModel(testModel);

        if sol.stat == 1
            growthMatrix(i,j) = sol.f;
        else
            growthMatrix(i,j) = NaN;
        end

    end

end


%% 8. Display sensitivity matrix

fprintf('\n============================================\n');
fprintf('Acetate / CO2 sensitivity matrix\n');
fprintf('============================================\n');

fprintf('Rows    = acetate uptake lower bounds\n');
fprintf('Columns = CO2 uptake lower bounds\n\n');

disp('CO2 levels:');
disp(co2Levels);

disp('Acetate levels:');
disp(acetateLevels);

disp('Growth matrix:');
disp(growthMatrix);


%% 9. Growth increase relative to baseline

growthIncrease = growthMatrix - baselineGrowth;

fprintf('\n============================================\n');
fprintf('Growth increase relative to baseline\n');
fprintf('============================================\n');

disp(growthIncrease);


%% 10. Find maximum biomass

[maxGrowth, linearIdx] = max(growthMatrix(:));

[rowIdx, colIdx] = ind2sub(size(growthMatrix), linearIdx);

fprintf('\n============================================\n');
fprintf('Maximum biomass production\n');
fprintf('============================================\n');

fprintf('Maximum biomass flux = %.6f\n', maxGrowth);
fprintf('Acetate LB           = %.2f\n', ...
    acetateLevels(rowIdx));
fprintf('CO2 LB               = %.2f\n', ...
    co2Levels(colIdx));


%% 11. Save results

resultsTable = array2table(growthMatrix);

% Add readable variable names
for j = 1:length(co2Levels)
    resultsTable.Properties.VariableNames{j} = ...
        sprintf('CO2_%g', co2Levels(j));
end

% Add acetate levels as first column
resultsTable = addvars(resultsTable, acetateLevels(:), ...
    'Before', 1, ...
    'NewVariableNames', 'Acetate_LB');

writetable(resultsTable, ...
    'Acetate_CO2_growth_sensitivity.xlsx');

fprintf('\nResults saved to:\n');
fprintf('Acetate_CO2_growth_sensitivity.xlsx\n');

