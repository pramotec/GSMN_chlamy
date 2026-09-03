%% ================= DIAMOND in MATLAB 自动化全流程 (完美兼容版) ================= %%
clear blast4414_Template bbh_data forward_data reverse_data;

% 1. 基础设置
query_fasta = '4414.fasta';
subject_fasta = 'iCre1355.fasta';
query_id = '4414';
subject_id = 'iCre1355';

% 2. 在 MATLAB 后台调用 DIAMOND 建库
disp('1/5 正在建立 DIAMOND 数据库...');
system(sprintf('diamond.exe makedb --in %s -d %s_db --quiet', subject_fasta, subject_id));
system(sprintf('diamond.exe makedb --in %s -d %s_db --quiet', query_fasta, query_id));

% 3. 在 MATLAB 后台运行 DIAMOND 双向比对
disp('2/5 正在运行 DIAMOND 双向高速比对 (请稍候)...');
system(sprintf('diamond.exe blastp -q %s -d %s_db -o forward.tsv --outfmt 6 --evalue 1e-10 --max-target-seqs 1 --quiet', query_fasta, subject_id));
system(sprintf('diamond.exe blastp -q %s -d %s_db -o reverse.tsv --outfmt 6 --evalue 1e-10 --max-target-seqs 1 --quiet', subject_fasta, query_id));

% 4. 读取 TSV 文件
disp('3/5 读取比对结果...');
varNames = {'query_id', 'subject_id', 'identity', 'alignment_length', ...
            'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'};
% 强制让 MATLAB 尝试自动识别数值类型，避免全部读成文本
opts = delimitedTextImportOptions('VariableNames', varNames, 'Delimiter', '\t', 'DataLines', 1);
forward_data = readtable('forward.tsv', opts);
reverse_data = readtable('reverse.tsv', opts);

% 5. 计算双向最佳比对 (BBH)
disp('4/5 正在计算双向最佳比对 (BBH)...');
[~, idx_f, ~] = unique(forward_data.query_id, 'stable');
f_best = forward_data(idx_f, :);

[~, idx_r, ~] = unique(reverse_data.query_id, 'stable');
r_best = reverse_data(idx_r, :);

map_reverse = containers.Map(r_best.query_id, r_best.subject_id);

bbh_idx = [];
for i = 1:height(f_best)
    q = f_best.query_id{i};
    s = f_best.subject_id{i};
    if isKey(map_reverse, s) && strcmp(map_reverse(s), q)
        bbh_idx(end+1) = i;
    end
end
bbh_data = f_best(bbh_idx, :);
fprintf('      => 共找到 %d 对双向最佳比对！\n', height(bbh_data));

% 6. 组装 RAVEN 所需的 blastStructure (包含所有数据类型和字段名补丁)
disp('5/5 正在封装 RAVEN 专属的 blastStructure...');
blast4414_Template = struct();

% 宏观 ID 补丁
blast4414_Template(1).fromId    = query_id;       
blast4414_Template(1).toId      = subject_id;   

% 微观基因补丁 (确保是 cellstr)
blast4414_Template(1).fromGenes = cellstr(bbh_data.query_id);   
blast4414_Template(1).toGenes   = cellstr(bbh_data.subject_id); 

% 数值类型强转补丁 (解决 '<' 报错)
% 使用 str2double(string(...)) 确保无论是文本还是数字都能转为 double 矩阵
blast4414_Template(1).evalue    = str2double(string(bbh_data.evalue));
blast4414_Template(1).identity  = str2double(string(bbh_data.identity));
blast4414_Template(1).aligLen   = str2double(string(bbh_data.alignment_length));
blast4414_Template(1).mismatch  = str2double(string(bbh_data.mismatches));

% 大小写补丁 (bitscore)
blast4414_Template(1).bitscore  = str2double(string(bbh_data.bit_score)); 

% 缺失字段补丁 (ppos)
blast4414_Template(1).ppos      = str2double(string(bbh_data.identity)); 

disp('==================================================');
disp('✅ 自动化流程执行完毕！');
disp('现在请执行以下命令开始重构：');
fprintf('TemplateCobra.id = ''%s'';\n', subject_id);
disp('Draft4414FromTemplate = getModelFromHomology({TemplateCobra}, blast4414_Template, ''4414'', {}, 2, false, 10^-5, 100);');