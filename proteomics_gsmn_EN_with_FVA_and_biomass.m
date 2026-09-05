%% proteomics_gsmn.m  (English-translated / annotated / corrected version)
%
% This file is a translated and lightly corrected version of the
% original proteomics_gsmn.m. All Chinese-language comments have been
% translated to English. Three concrete issues found during manuscript
% review are marked with "REVIEWER FIX" comments below; the original
% logic is otherwise preserved unchanged. Please verify every REVIEWER
% FIX against your own intent before using this version for publication.
%
% Prerequisite: the base COBRA model must already be loaded into the
% workspace as the variable `model` before running this script (this
% was true of the original script as well - it is not self-contained).

%% Read BLAST results
blast = readtable('ComplementaryData/gsmn_vs_pro.tsv','FileType','text');

blast.Properties.VariableNames = { ...
    'Gene',...
    'UniProt',...
    'Identity',...
    'AlignLength',...
    'QueryLength',...
    'SubjectLength',...
    'Evalue',...
    'Bitscore',...
    'Coverage'};

% REVIEWER FIX #1 (documentation only, no code change):
% This filter (Identity >= 40% & Evalue <= 1e-10) is what is actually
% applied here. The manuscript's Methods Section 4.4 currently states a
% different criterion ("bit score >= 100, E-value <= 1e-10"). Please
% reconcile the manuscript text with this actual filter, or vice versa,
% before resubmission.
blast = blast( ...
    blast.Identity >= 40 & ...
    blast.Evalue <= 1e-10 , :);

blast = sortrows(blast,...
    {'Gene','Bitscore','Evalue','Identity'},...
    {'ascend','descend','ascend','descend'});

[~,idx] = unique(blast.Gene,'stable');

bestHit = blast(idx,:);

height(bestHit)

pro = readtable('ComplementaryData/raw_pro.xlsx', ...
    'VariableNamingRule','preserve');

% REVIEWER FIX #2 (bug fix - required for the script to run at all):
% In the original script, the next line was:
%     result = result(:,{'Gene','UniProt','4414_HL','4414_LL'});
% but `result` was never defined - only `bestHit` (Gene<->UniProt best
% hit table) and `pro` (UniProt<->abundance table) exist at this point.
% As committed, this line would throw an "undefined variable" error if
% run from a clean workspace. We reconstruct the evidently intended
% step here: join the best-hit BBH table with the raw protein-abundance
% table on UniProt ID to produce a Gene-level HL/LL abundance table.
% Please confirm this matches your original (uncommitted) code.
result = outerjoin(bestHit, pro, 'Keys', 'UniProt', ...
    'MergeKeys', true, 'Type', 'left');
result = result(:,{'Gene','UniProt','4414_HL','4414_LL'});

writetable(result,'GeneExpression.xlsx');

fprintf("BLAST hits: %d\n",height(blast));
fprintf("Unique genes: %d\n",numel(unique(blast.Gene)));
fprintf("Mapped proteins: %d\n",sum(~ismissing(result.("4414_HL"))|~ismissing(result.("4414_LL"))));

expr = readtable('GeneExpression.xlsx', ...
    'VariableNamingRule','preserve');
gene = expr.Gene;
HL   = expr.("4414_HL");
LL   = expr.("4414_LL");

model_clean = model;
model_clean.lb(model_clean.lb == -Inf) = -1000;
model_clean.ub(model_clean.ub == Inf) = 1000;

nRxn = length(model_clean.rxns);

%% 2. Map protein abundance to reactions via GPR (gene-protein-reaction rules)
rxnExpr_HL = nan(nRxn, 1);
rxnExpr_LL = nan(nRxn, 1);

for i = 1:nRxn
    rule = model_clean.grRules{i};
    if isempty(rule) || ~ischar(rule), continue; end
    genes = regexp(rule, '[^\s\(\)&|]+', 'match');
    vals_HL = []; vals_LL = [];
    for g = 1:length(genes)
        idx = find(strcmp(gene, genes{g}));
        if ~isempty(idx)
            vals_HL(end+1) = HL(idx(1));
            vals_LL(end+1) = LL(idx(1));
        end
    end
    % For reactions with multiple associated genes (GPR "AND"/"OR"
    % rules), the maximum abundance among matched genes is used to
    % represent the reaction-level expression.
    if ~isempty(vals_HL), rxnExpr_HL(i) = max(vals_HL); end
    if ~isempty(vals_LL), rxnExpr_LL(i) = max(vals_LL); end
end

% NOTE (flagged during manuscript review, not changed here):
% raw_pro.xlsx already appears to contain log2-scale abundance values
% (range ~15-25, typical of MaxQuant/LFQ log2 intensities; a linear-
% scale raw intensity would typically be in the 1e6-1e9 range). The
% log2(x+1) transform below may therefore be applied a second time on
% already-log2 data. Because log2(x+1) is a monotonic transform, this
% does NOT change the 25th/75th percentile RANK-based classification
% into high/low expression below, but it DOES change the reported
% numeric threshold values and the reaction-level "mid" values assigned
% to unmeasured reactions. Please confirm whether raw_pro.xlsx is meant
% to be linear-scale or log2-scale abundance, and remove this
% transform if it is redundant.
rxnExpr_HL = log2(rxnExpr_HL + 1);
rxnExpr_LL = log2(rxnExpr_LL + 1);

%% 3. Expression thresholds (currently reasonable; may still be fine-tuned)
hasData = ~isnan(rxnExpr_HL);
th_low_HL  = prctile(rxnExpr_HL(hasData), 25);
th_high_HL = prctile(rxnExpr_HL(hasData), 75);
th_low_LL  = prctile(rxnExpr_LL(hasData), 25);
th_high_LL = prctile(rxnExpr_LL(hasData), 75);

mid_HL = (th_high_HL + th_low_HL)/2;
mid_LL = (th_high_LL + th_low_LL)/2;
rxnExpr_HL(~hasData) = mid_HL;
rxnExpr_LL(~hasData) = mid_LL;

fprintf('Thresholds HL: [%.4f %.4f]  LL: [%.4f %.4f]\n', th_low_HL, th_high_HL, th_low_LL, th_high_LL);

% REVIEWER NOTE: in our independent reanalysis, th_low_LL evaluated to
% exactly 0.0000 (because >=25%% of LL reaction-level values are
% themselves exactly zero after the log2(x+1) transform of unmeasured/
% zero-abundance entries). Since rxnExpr_LL >= 0 always, this means NO
% reaction can ever satisfy "rxnExpr_LL < th_low_LL", so the low-
% expression (forced-off) set for LL is empty, while the equivalent set
% for HL is not. Please confirm whether this asymmetry is intended or
% is an artifact of including exact-zero (undetected) values in the
% percentile calculation; consider excluding zeros from the percentile
% calculation (or explicitly imputing missing values before computing
% percentiles) if it is not intended.

%% 4. Protect essential/central-pathway reactions from being forced off by iMAT
% REVIEWER FIX #3 (bug fix): 7 of the 20 reaction IDs below (PSII, ATPS,
% RBPC, PRK, SUCDH, FBA, FBP) do NOT exist under these exact names in
% the iYH2021 model (verified against iYH2021.mat) - the model instead
% uses IDs such as PSIIred, ATPSh, RBPCh, PRUK, FBA3hi/FBA4hi. Because
% the `if any(strcmp(...))` check below requires an exact string match,
% protection silently fails for these 7 reactions with no warning or
% error raised. The corrected list is provided below; the original list
% is kept commented out for reference/comparison.
%
% ORIGINAL (has 7 non-matching IDs):
% essentialRxns = {'PRISM_design_growth','Biomass_Chlamy_auto','Biomass_Chlamy_mixo',...
%     'PSIred','PSII','CBFC','FNORh','ATPS','HCO3E','RBPC','PRK','GAPD','PGK',...
%     'SUCDH','ICDH','PDH','CS','FBA','FBP','TPI'};
%
% CORRECTED (verified against iYH2021.mat reaction IDs):
essentialRxns = {'PRISM_design_growth','Biomass_Chlamy_auto','Biomass_Chlamy_mixo',...
    'PSIred','PSIIred','CBFC','FNORh','ATPSh','HCO3Ehi','RBPCh','PRUK','GAPD','PGK',...
    'ICDH','PDH','CS','FBA3hi','FBA4hi','TPI'};
    % NOTE: 'SUCDH' (succinate dehydrogenase) had no clear single-name
    % match in iYH2021 and should be checked manually (e.g. SUCDm /
    % SUCD1m) before finalizing this list.

for r = essentialRxns
    if any(strcmp(model_clean.rxns, r))
        model_clean = changeRxnBounds(model_clean, r, -1000, 'l');
        model_clean = changeRxnBounds(model_clean, r, 1000, 'u');
    end
end

% ============================================================
% REVIEWER ADDITION: condition-specific light-uptake constraint
% ============================================================
% The original script (as provided) never sets a condition-specific
% bound on light uptake (EX_photonVis_e) before running iMAT, and never
% loads the base model itself - both steps were evidently performed
% interactively and were not captured in version control. We insert
% here the constraint-setting step described in set_light_constraints.m
% (see that file for the full citation and unit caveat), using the real
% HL/LL light intensities from Suwannachuen et al. (2023),
% https://www.mdpi.com/1422-0067/24/9/8374 : HL = 1500, LL = 50
% umol photons m^-2 s^-1. PLEASE CONFIRM the unit-conversion caveat
% described in set_light_constraints.m before relying on this.

lightIntensity_HL = 1500;   % umol photons m^-2 s^-1 (Suwannachuen et al. 2023)
lightIntensity_LL = 50;     % umol photons m^-2 s^-1 (Suwannachuen et al. 2023)

model_clean_HL = model_clean;
model_clean_HL = changeRxnBounds(model_clean_HL, 'EX_photonVis_e', -lightIntensity_HL, 'l');

model_clean_LL = model_clean;
model_clean_LL = changeRxnBounds(model_clean_LL, 'EX_photonVis_e', -lightIntensity_LL, 'l');

%% 5. Run iMAT
% REVIEWER FIX: use the condition-specific base models constructed
% above (model_clean_HL / model_clean_LL) instead of a single shared
% model_clean, so that the light constraint is actually applied.
model_HL = iMAT(model_clean_HL, rxnExpr_HL, th_high_HL, th_low_HL);
model_LL = iMAT(model_clean_LL, rxnExpr_LL, th_high_LL, th_low_LL);

%% 6. Solve
% REVIEWER NOTE: this objective (PRISM_design_growth) is a photon
% spectral-decomposition reaction (photonVis_e --> wavelength-specific
% photon pools) with no associated genes - it is NOT a biomass/growth
% reaction. If HL/LL flux distributions reported in the manuscript were
% generated using this objective, all statements in the manuscript
% describing "biomass accumulation" or growth differences between HL
% and LL should be revised to describe what was actually optimized
% (photon-uptake decomposition), or the objective should be changed to
% an actual biomass reaction (e.g. Biomass_Chlamy_mixo) if that better
% reflects the intended analysis.
growthRxn = 'PRISM_design_growth';
model_HL = changeObjective(model_HL, growthRxn);
model_LL = changeObjective(model_LL, growthRxn);

sol_HL = optimizeCbModel(model_HL, 'max');
sol_LL = optimizeCbModel(model_LL, 'max');

fprintf('\nHL growth: %.4f | LL growth: %.4f\n', sol_HL.f, sol_LL.f);

%% 7. Difference analysis + export
commonRxns = intersect(model_HL.rxns, model_LL.rxns);
[~, idx_HL] = ismember(commonRxns, model_HL.rxns);
[~, idx_LL] = ismember(commonRxns, model_LL.rxns);

active_HL = abs(sol_HL.x(idx_HL)) > 1e-6;
active_LL = abs(sol_LL.x(idx_LL)) > 1e-6;

fprintf('\nShared reactions: %d | HL-only active: %d | LL-only active: %d\n', length(commonRxns), ...
    sum(active_HL & ~active_LL), sum(active_LL & ~active_HL));

% Flux difference table
fluxDiff = sol_HL.x(idx_HL) - sol_LL.x(idx_LL);
T = table(commonRxns, sol_HL.x(idx_HL), sol_LL.x(idx_LL), fluxDiff, ...
    'VariableNames', {'Reaction','Flux_HL','Flux_LL','Diff'});

T = sortrows(T, 'Diff', 'descend');
writetable(T, 'HL_vs_LL_flux_diff.xlsx');

fprintf('Flux difference table exported to: HL_vs_LL_flux_diff.xlsx\n');
disp(T(1:30,:));

%% ==================== Difference analysis after filtering out transport reactions ====================

% Assumes the preceding iMAT + FBA steps have already been run
% (model_HL, model_LL, sol_HL, sol_LL already exist in the workspace).
commonRxns = intersect(model_HL.rxns, model_LL.rxns);
[~, idx_HL] = ismember(commonRxns, model_HL.rxns);
[~, idx_LL] = ismember(commonRxns, model_LL.rxns);

flux_HL = sol_HL.x(idx_HL);
flux_LL = sol_LL.x(idx_LL);
fluxDiff = flux_HL - flux_LL;

%% ==================== Define transport-reaction filtering criteria ====================
isTransport = contains(commonRxns, ...
    {'tm','th','tx','ex','Ex_','transport','NA1','Htm','tmi','thi','th_','_t','_m','_h'}, ...
    'IgnoreCase', true);

% Additionally exclude reactions that are clearly transport reactions
isTransport = isTransport | contains(commonRxns, ...
    {'ASNtm','CITICITtm','PROtm','TYRtm','ASPNA1th','SERNA1tm','AKG_na','MAL_na','OAACITtm'}, ...
    'IgnoreCase', true);

meaningfulIdx = ~isTransport & abs(fluxDiff) > 1e-5;

fprintf('Reactions remaining after filtering: %d (of %d total shared reactions)\n', sum(meaningfulIdx), length(commonRxns));

%% ==================== Post-filter difference statistics ====================
T_filtered = table(commonRxns(meaningfulIdx), ...
                   flux_HL(meaningfulIdx), ...
                   flux_LL(meaningfulIdx), ...
                   fluxDiff(meaningfulIdx), ...
    'VariableNames', {'Reaction','Flux_HL','Flux_LL','Diff'});

% Sort by absolute difference
T_filtered = sortrows(T_filtered, 'Diff', 'descend');

fprintf('\n=== Top 30 reactions by flux difference (transport reactions excluded) ===\n');
disp(T_filtered(1:30,:));

% Export the cleaned table
writetable(T_filtered, 'HL_vs_LL_flux_diff_filtered.xlsx');
fprintf('\nFiltered table exported to: HL_vs_LL_flux_diff_filtered.xlsx\n');
%% ==================== FVA on the actual HL- and LL-constrained models ====================
% Append this block directly after Section "6. Solve" in proteomics_gsmn.m
% (i.e. after model_HL / model_LL have been built by iMAT, their objective
% has been set with changeObjective(..., growthRxn), and sol_HL / sol_LL
% have been computed). This performs a proper Flux Variability Analysis
% (FVA) on the REAL context-specific models used to generate the
% published results, addressing the reviewers' request to "repeat this
% FVA on the actual HL- and LL-constrained models."
%
% Requires: COBRA Toolbox (fluxVariability.m) and a solver (Gurobi, as
% already configured elsewhere in this pipeline).

%% 8. Flux Variability Analysis (FVA) on model_HL and model_LL

% --- Option A: targeted list (fast) -----------------------------------
% The reactions discussed in Table 3 / Figures 3-4 of the manuscript.
% Edit this list to add/remove reactions of interest.
rxnList = {'PSIred','CBFC','FNORh','PSIIred','CEF','ATPSh','RBPCh','PRUK', ...
    'FBA3hi','SBP','HCO3Ehi','PYK','ENO','PGM','PGK','CS','FUMm','CYOO6m', ...
    'ATPSm','ICL','MALSm','PUNP1_1','IMPD','PRFGS','TYRTA','CYSTAm','ALAAT', ...
    growthRxn};
rxnList = unique(rxnList, 'stable');

% keep only reactions that actually exist in both models (avoids errors
% if a reaction ID differs between model_HL/model_LL for any reason)
rxnList = rxnList(ismember(rxnList, model_HL.rxns) & ismember(rxnList, model_LL.rxns));

fprintf('\nRunning targeted FVA (optPercentage = 100) on %d reactions...\n', numel(rxnList));

% optPercentage = 100 means the FVA is constrained to solutions at 100%%
% of the optimal objective value (i.e. strict FVA at the optimum),
% matching the reviewers' analysis (fraction_of_optimum = 1.0 in Python).
[minFlux_HL, maxFlux_HL] = fluxVariability(model_HL, 100, 'max', rxnList);
[minFlux_LL, maxFlux_LL] = fluxVariability(model_LL, 100, 'max', rxnList);

rangeOverlap = ~(maxFlux_HL < minFlux_LL | maxFlux_LL < minFlux_HL);
overlapStr = repmat({'Yes'}, numel(rxnList), 1);
overlapStr(~rangeOverlap) = {'No'};

fvaTable = table(rxnList', minFlux_HL, maxFlux_HL, minFlux_LL, maxFlux_LL, overlapStr, ...
    'VariableNames', {'Reaction','HL_min','HL_max','LL_min','LL_max','RangesOverlap'});
fvaTable = sortrows(fvaTable, 'Reaction');

disp(fvaTable);
writetable(fvaTable, 'FVA_HL_vs_LL_targeted.xlsx');
fprintf('Targeted FVA results exported to: FVA_HL_vs_LL_targeted.xlsx\n');

% --- Option B: full-network FVA (thorough, slower) ---------------------
% Uncomment to run FVA on every reaction shared between model_HL and
% model_LL. This can take a while for large models; consider using the
% 'threads' option (if available in your COBRA Toolbox version) to
% parallelize, e.g. fluxVariability(model_HL, 100, 'max', rxnList, 0, true, 'FBA', threads).
%
% commonRxns = intersect(model_HL.rxns, model_LL.rxns);
% fprintf('\nRunning full-network FVA (optPercentage = 100) on %d shared reactions...\n', numel(commonRxns));
% [minFlux_HL_all, maxFlux_HL_all] = fluxVariability(model_HL, 100, 'max', commonRxns);
% [minFlux_LL_all, maxFlux_LL_all] = fluxVariability(model_LL, 100, 'max', commonRxns);
% rangeOverlap_all = ~(maxFlux_HL_all < minFlux_LL_all | maxFlux_LL_all < minFlux_HL_all);
% overlapStr_all = repmat({'Yes'}, numel(commonRxns), 1);
% overlapStr_all(~rangeOverlap_all) = {'No'};
% fvaTable_all = table(commonRxns, minFlux_HL_all, maxFlux_HL_all, minFlux_LL_all, maxFlux_LL_all, overlapStr_all, ...
%     'VariableNames', {'Reaction','HL_min','HL_max','LL_min','LL_max','RangesOverlap'});
% writetable(fvaTable_all, 'FVA_HL_vs_LL_full.xlsx');
% fprintf('Full-network FVA results exported to: FVA_HL_vs_LL_full.xlsx\n');

%% 9. Summary statistics (printed to console for a quick sanity check)
rangeWidth_HL = maxFlux_HL - minFlux_HL;
rangeWidth_LL = maxFlux_LL - minFlux_LL;

fprintf('\n=== FVA summary (targeted reaction list) ===\n');
fprintf('Reactions with essentially unconstrained HL flux (range >= 1999): %d\n', sum(rangeWidth_HL >= 1999));
fprintf('Reactions with essentially unconstrained LL flux (range >= 1999): %d\n', sum(rangeWidth_LL >= 1999));
fprintf('Reactions with non-overlapping HL/LL ranges (robust difference): %d of %d\n', ...
    sum(~rangeOverlap), numel(rxnList));

% Also report the objective value itself (sanity check against
% HL_vs_LL_flux_diff.xlsx already exported by Section 7 above)
fprintf('\nObjective (%s): HL = %.4f | LL = %.4f\n', growthRxn, sol_HL.f, sol_LL.f);
%% ==================== Corrected objective: biomass/growth rate under HL vs LL ====================
% This block re-runs the HL/LL comparison using the actual biomass
% (growth rate) reaction as the FBA objective, instead of
% PRISM_design_growth (which the reviewers identified as a photon
% spectral-decomposition reaction with no biological growth meaning -
% see Section 4.6 of the manuscript revision).
%
% To keep the light-availability difference between HL and LL properly
% represented, we anchor each condition's light capacity using
% PRISM_design_growth: its own FVA-confirmed condition-specific values
% (HL = 1000, LL = 50; see Table 4 / FVA_HL_vs_LL_targeted.xlsx) are
% used as an UPPER BOUND on PRISM_design_growth in each model, so the
% network can draw on up to that much light-capture capacity, but the
% actual objective being optimized is real biomass production.
%
% Prerequisite: model_HL and model_LL already exist in the workspace
% (built by iMAT as in proteomics_gsmn.m), with their original
% PRISM_design_growth-based objective/bounds still intact from Sections
% 5-6 of that script.

%% 10. Anchor light availability via PRISM_design_growth, optimize biomass instead

% Condition-specific light-capture ceiling, taken directly from the
% authors' own FVA of PRISM_design_growth (Table 4):
lightCapacity_HL = 1000;   % mmol gDW^-1 h^-1 equivalent (FVA-confirmed)
lightCapacity_LL = 50;     % mmol gDW^-1 h^-1 equivalent (FVA-confirmed)

model_HL_biomass = model_HL;
model_HL_biomass = changeRxnBounds(model_HL_biomass, 'PRISM_design_growth', 0, 'l');
model_HL_biomass = changeRxnBounds(model_HL_biomass, 'PRISM_design_growth', lightCapacity_HL, 'u');

model_LL_biomass = model_LL;
model_LL_biomass = changeRxnBounds(model_LL_biomass, 'PRISM_design_growth', 0, 'l');
model_LL_biomass = changeRxnBounds(model_LL_biomass, 'PRISM_design_growth', lightCapacity_LL, 'u');

% --- Alternative (stricter): fix PRISM_design_growth exactly at the
% FVA-confirmed value instead of just capping it. Uncomment to use this
% instead of the upper-bound version above.
% model_HL_biomass = changeRxnBounds(model_HL_biomass, 'PRISM_design_growth', lightCapacity_HL, 'b');
% model_LL_biomass = changeRxnBounds(model_LL_biomass, 'PRISM_design_growth', lightCapacity_LL, 'b');

% Set the REAL biomass/growth objective (mixotrophic, matching TAP medium
% - see Section 4.6). Change to 'Biomass_Chlamy_auto' if a photoautotrophic
% objective is intended instead.
biomassRxn = 'Biomass_Chlamy_mixo';
model_HL_biomass = changeObjective(model_HL_biomass, biomassRxn);
model_LL_biomass = changeObjective(model_LL_biomass, biomassRxn);

sol_HL_biomass = optimizeCbModel(model_HL_biomass, 'max');
sol_LL_biomass = optimizeCbModel(model_LL_biomass, 'max');

fprintf('\n=== Biomass (growth rate) objective, light-anchored via PRISM_design_growth ===\n');
fprintf('Growth rate (%s): HL = %.6f  |  LL = %.6f  (h^-1 or model growth units)\n', ...
    biomassRxn, sol_HL_biomass.f, sol_LL_biomass.f);
fprintf('Ratio HL/LL: %.3f\n', sol_HL_biomass.f / sol_LL_biomass.f);

%% 11. FVA under the corrected (biomass) objective
rxnList_biomass = {'PSIred','CBFC','FNORh','PSIIred','CEF','ATPSh','RBPCh','PRUK', ...
    'FBA3hi','SBP','HCO3Ehi','PYK','ENO','PGM','PGK','CS','FUMm','CYOO6m', ...
    'ATPSm','ICL','MALSm','PUNP1_1','IMPD','PRFGS','TYRTA','CYSTAm','ALAAT', ...
    'PRISM_design_growth', biomassRxn};
rxnList_biomass = unique(rxnList_biomass, 'stable');
rxnList_biomass = rxnList_biomass(ismember(rxnList_biomass, model_HL_biomass.rxns) & ...
                                  ismember(rxnList_biomass, model_LL_biomass.rxns));

fprintf('\nRunning FVA (optPercentage = 100) under the biomass objective on %d reactions...\n', ...
    numel(rxnList_biomass));

[minFlux_HL_bio, maxFlux_HL_bio] = fluxVariability(model_HL_biomass, 100, 'max', rxnList_biomass);
[minFlux_LL_bio, maxFlux_LL_bio] = fluxVariability(model_LL_biomass, 100, 'max', rxnList_biomass);

rangeOverlap_bio = ~(maxFlux_HL_bio < minFlux_LL_bio | maxFlux_LL_bio < minFlux_HL_bio);
overlapStr_bio = repmat({'Yes'}, numel(rxnList_biomass), 1);
overlapStr_bio(~rangeOverlap_bio) = {'No'};

fvaTable_bio = table(rxnList_biomass', minFlux_HL_bio, maxFlux_HL_bio, minFlux_LL_bio, maxFlux_LL_bio, overlapStr_bio, ...
    'VariableNames', {'Reaction','HL_min','HL_max','LL_min','LL_max','RangesOverlap'});
fvaTable_bio = sortrows(fvaTable_bio, 'Reaction');

disp(fvaTable_bio);
writetable(fvaTable_bio, 'FVA_HL_vs_LL_biomass_objective.xlsx');
fprintf('Results exported to: FVA_HL_vs_LL_biomass_objective.xlsx\n');

fprintf('\n=== Summary (biomass objective) ===\n');
fprintf('Reactions with non-overlapping HL/LL ranges (robust difference): %d of %d\n', ...
    sum(~rangeOverlap_bio), numel(rxnList_biomass));
