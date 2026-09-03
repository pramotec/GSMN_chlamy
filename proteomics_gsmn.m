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
result = result(:,{'Gene','UniProt','4414_HL','4414_LL'});

writetable(result,'GeneExpression.xlsx');

fprintf("BLAST hits: %d\n",height(blast));
fprintf("Unique genes: %d\n",numel(unique(blast.Gene)));
fprintf("Mapped proteins: %d\n",sum(~ismissing(result.4414_HL)|~ismissing(result.4414_LL)));

expr = readtable('GeneExpression.xlsx', ...
    'VariableNamingRule','preserve');


model_clean = model;
model_clean.lb(model_clean.lb == -Inf) = -1000;
model_clean.ub(model_clean.ub == Inf) = 1000;

nRxn = length(model_clean.rxns);

%% 2. 蛋白质组处理
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
    if ~isempty(vals_HL), rxnExpr_HL(i) = max(vals_HL); end
    if ~isempty(vals_LL), rxnExpr_LL(i) = max(vals_LL); end
end

rxnExpr_HL = log2(rxnExpr_HL + 1);
rxnExpr_LL = log2(rxnExpr_LL + 1);

%% 3. 阈值（当前较好，可继续微调）
hasData = ~isnan(rxnExpr_HL);
th_low_HL  = prctile(rxnExpr_HL(hasData), 25);
th_high_HL = prctile(rxnExpr_HL(hasData), 75);
th_low_LL  = prctile(rxnExpr_LL(hasData), 25);
th_high_LL = prctile(rxnExpr_LL(hasData), 75);

mid_HL = (th_high_HL + th_low_HL)/2;
mid_LL = (th_high_LL + th_low_LL)/2;
rxnExpr_HL(~hasData) = mid_HL;
rxnExpr_LL(~hasData) = mid_LL;

fprintf('阈值 HL: [%.4f %.4f] LL: [%.4f %.4f]\n', th_low_HL, th_high_HL, th_low_LL, th_high_LL);

%% 4. 保护关键反应
essentialRxns = {'PRISM_design_growth','Biomass_Chlamy_auto','Biomass_Chlamy_mixo',...
    'PSIred','PSII','CBFC','FNORh','ATPS','HCO3E','RBPC','PRK','GAPD','PGK',...
    'SUCDH','ICDH','PDH','CS','FBA','FBP','TPI'};

for r = essentialRxns
    if any(strcmp(model_clean.rxns, r))
        model_clean = changeRxnBounds(model_clean, r, -1000, 'l');
        model_clean = changeRxnBounds(model_clean, r, 1000, 'u');
    end
end

%% 5. 运行 iMAT
model_HL = iMAT(model_clean, rxnExpr_HL, th_high_HL, th_low_HL);
model_LL = iMAT(model_clean, rxnExpr_LL, th_high_LL, th_low_LL);

%% 6. 求解
growthRxn = 'PRISM_design_growth';
model_HL = changeObjective(model_HL, growthRxn);
model_LL = changeObjective(model_LL, growthRxn);

sol_HL = optimizeCbModel(model_HL, 'max');
sol_LL = optimizeCbModel(model_LL, 'max');

fprintf('\nHL growth: %.4f | LL growth: %.4f\n', sol_HL.f, sol_LL.f);

%% 7. 差异分析 + 导出
commonRxns = intersect(model_HL.rxns, model_LL.rxns);
[~, idx_HL] = ismember(commonRxns, model_HL.rxns);
[~, idx_LL] = ismember(commonRxns, model_LL.rxns);

active_HL = abs(sol_HL.x(idx_HL)) > 1e-6;
active_LL = abs(sol_LL.x(idx_LL)) > 1e-6;

fprintf('\n共有反应: %d | 只HL: %d | 只LL: %d\n', length(commonRxns), ...
    sum(active_HL & ~active_LL), sum(active_LL & ~active_HL));

% 通量差异表
fluxDiff = sol_HL.x(idx_HL) - sol_LL.x(idx_LL);
T = table(commonRxns, sol_HL.x(idx_HL), sol_LL.x(idx_LL), fluxDiff, ...
    'VariableNames', {'Reaction','Flux_HL','Flux_LL','Diff'});

T = sortrows(T, 'Diff', 'descend');
writetable(T, 'HL_vs_LL_flux_diff.xlsx');

fprintf('差异表格已导出到: HL_vs_LL_flux_diff.xlsx\n');
disp(T(1:30,:));   

%% ==================== 过滤转运反应后的差异分析 ====================

% 假设你已经运行完前面的 iMAT + FBA（model_HL, model_LL, sol_HL, sol_LL 已存在）
commonRxns = intersect(model_HL.rxns, model_LL.rxns);
[~, idx_HL] = ismember(commonRxns, model_HL.rxns);
[~, idx_LL] = ismember(commonRxns, model_LL.rxns);

flux_HL = sol_HL.x(idx_HL);
flux_LL = sol_LL.x(idx_LL);
fluxDiff = flux_HL - flux_LL;

%% ==================== 定义转运反应过滤条件 ====================
isTransport = contains(commonRxns, ...
    {'tm','th','tx','ex','Ex_','transport','NA1','Htm','tmi','thi','th_','_t','_m','_h'}, ...
    'IgnoreCase', true);

% 进一步排除明显转运
isTransport = isTransport | contains(commonRxns, ...
    {'ASNtm','CITICITtm','PROtm','TYRtm','ASPNA1th','SERNA1tm','AKG_na','MAL_na','OAACITtm'}, ...
    'IgnoreCase', true);

meaningfulIdx = ~isTransport & abs(fluxDiff) > 1e-5;

fprintf('过滤后剩余有意义反应数: %d (原共有 %d)\n', sum(meaningfulIdx), length(commonRxns));

%% ==================== 过滤后差异统计 ====================
T_filtered = table(commonRxns(meaningfulIdx), ...
                   flux_HL(meaningfulIdx), ...
                   flux_LL(meaningfulIdx), ...
                   fluxDiff(meaningfulIdx), ...
    'VariableNames', {'Reaction','Flux_HL','Flux_LL','Diff'});

% 按差异绝对值排序
T_filtered = sortrows(T_filtered, 'Diff', 'descend');

fprintf('\n=== 过滤转运后 —— 通量差异最大的前 30 个反应 ===\n');
disp(T_filtered(1:30,:));

% 导出干净表格
writetable(T_filtered, 'HL_vs_LL_flux_diff_filtered.xlsx');
fprintf('\n已导出过滤版表格: HL_vs_LL_flux_diff_filtered.xlsx\n');