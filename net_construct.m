%% STEP 1: TEMPLATE MODEL PREPARATION
% use iCre1355 as template model
initCobraToolbox;
iCreCobra = readCbModel('ComplementaryData/iCre1355_mixo.xml'); 
iCreRaven = ravenCobraWrapper(iCreCobra);

for i = 1:numel(iCreRaven.metComps)
    if iCreRaven.metComps(i) == 0;
        iCreRaven.metComps(i) = 1;
    else
    end
end
%% STEP 2: SEQUENCE ANALYSIS FOR CONSTRUCTING A METABOLIC NETWORK
iCreRaven.id = 'iCre';

blastAniCre=getBlast('An','ComplementaryData/4414.fasta','iCre','ComplementaryData/iCre1355.fasta');

save('ComplementaryData/blastAniCre.mat','blastAniCre');

load('ComplementaryData/blastAniCre.mat');

AnDraftFromiCre=getModelFromHomology({iCreRaven},blastAniCre,'An',{},2,false,10^-5,100);

rxnToAdd.rxns =setdiff(iCreRaven.rxns,AnDraftFromiCre.rxns);

% comps
oldComps_An = AnDraftFromiCre.comps;
newComps    = iCreRaven.comps;

[~, map_An] = ismember(oldComps_An, newComps);
AnDraftFromiCre.metComps = map_An(AnDraftFromiCre.metComps);
AnDraftFromiCre.comps = newComps;

if isfield(iCreRaven,'compNames')
    AnDraftFromiCre.compNames = iCreRaven.compNames;
end

AnDraftPlusiCre = addRxnsGenesMets(AnDraftFromiCre,iCreRaven,rxnToAdd.rxns,...
    false,'additional rxns from iCre1355 model',2);

goodAnDraftPlusiCre = AnDraftPlusiCre

goodAnDraftPlusiCre = setParam(goodAnDraftPlusiCre,'obj',{'Biomass_Chlamy_mixo'},1);
sol = solveLP(goodAnDraftPlusiCre,1);
printFluxes(goodAnDraftPlusiCre, sol.x,true);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);

model_An=contractModel(goodAnDraftPlusiCre);

printModelStats(model_An,true,true);
%exportForGit(model_An,'model_An','');

%% STEP 3: GAP-FILLING AND POLISHING GSMM 
% 3.1 update grRules
model_Anplus = model_An;
newGrRules = readcell('ComplementaryData/updated_4414.xlsx', 'Sheet', 'updated_grRules');
rxnID = newGrRules(2:end, 1);
geneAssoc = newGrRules(2:end, 2);
[ismem, rxnIdx] = ismember(rxnID, model_Anplus.rxns);

for i = 1:numel(rxnID)
    if ismem(i)
        model_Anplus.grRules{rxnIdx(i)} = geneAssoc{i};
    else
        warning(['no reaction', rxnID{i}]);
    end
end
model_Anplus = generateRules(model_Anplus);

% 3.1 update metabolite names
model_Anplus2=model_Anplus;
[~, textData1]=xlsread('ComplementaryData/updated_4414.xlsx','updated_MetNames');
metNames = struct();
metNames.old = textData1(2:end,1);
metNames.new = textData1(2:end,2);
[a, b]=ismember(model_Anplus2.metNames,metNames.old);
I=find(a);
model_Anplus2.metNames(I)=metNames.new(b(I));
printModelStats(model_Anplus2,true,true);

% 3.1 update EC numbers
model_Anplus3=model_Anplus2;
[~, textData1]=xlsread('ComplementaryData/updated_4414.xlsx','updated_RxnECs');
updated_RxnECs = struct();
updated_RxnECs.rxns = textData1(2:end,1);
updated_RxnECs.new = textData1(2:end,2);
[a, b]=ismember(model_Anplus3.rxns,updated_RxnECs.rxns);
I=find(a);
model_Anplus3.eccodes(I)=updated_RxnECs.new(b(I));
invalidIdx = cellfun(@(x) ~ischar(x) && ~isstring(x), model_Anplus3.eccodes);
model_Anplus3.eccodes(invalidIdx) = {''};

model_Anplus3.eccodes = model_Anplus3.eccodes(:);

numRxns = numel(model_Anplus3.rxns);
if numel(model_Anplus3.eccodes) < numRxns
    model_Anplus3.eccodes(end+1 : numRxns, 1) = {''};
end
model_Anplus3=contractModel(model_Anplus3);

printModelStats(model_Anplus3,true,true);

% 4 add new reactions from CreaveMe model

[~, SheetS] = xlsread('ComplementaryData/updated_4414.xlsx','new_Reactions2');
an_newRxns2 = struct();
an_newRxns2.rxns = SheetS(2:end,1);
an_newRxns2.rxnNames = SheetS(2:end,2);
an_newRxns2.equations = SheetS(2:end,3);
an_newRxns2.eccodes = SheetS(2:end,4);
an_newRxns2.grRules = SheetS(2:end,5);
an_newRxns2.subSystems = SheetS(2:end,6);
new1 = model_Anplus3
new2 = addRxns(new1, an_newRxns2, 3, '', true, true);
for i = 1:numel(an_newRxns2.rxns)
    new2 = setParam(new2, 'ub', an_newRxns2.rxns{i}, 1000);
end

model_final = new2;

genesToRemove = {
    'OLD_iCre_ChreCp001'; 'OLD_iCre_ChreCp002'; 'OLD_iCre_ChreCp003';
    'OLD_iCre_ChreCp004'; 'OLD_iCre_ChreCp008'; 'OLD_iCre_ChreCp009';
    'OLD_iCre_ChreCp019'; 'OLD_iCre_ChreCp021'; 'OLD_iCre_ChreCp023';
    'OLD_iCre_ChreCp026'; 'OLD_iCre_ChreCp027'; 'OLD_iCre_ChreCp029';
    'OLD_iCre_ChreCp031'; 'OLD_iCre_ChreCp032'; 'OLD_iCre_ChreCp040';
    'OLD_iCre_ChreCp043'; 'OLD_iCre_ChreCp044'; 'OLD_iCre_ChreCp045';
    'OLD_iCre_ChreCp048'; 'OLD_iCre_ChreCp049'; 'OLD_iCre_ChreCp050';
    'OLD_iCre_ChreCp051'; 'OLD_iCre_ChreCp053'; 'OLD_iCre_ChreCp054';
    'OLD_iCre_ChreCp056'; 'OLD_iCre_ChreCp057'; 'OLD_iCre_ChreCp058';
    'OLD_iCre_ChreCp061'; 'OLD_iCre_ChreCp062'; 'OLD_iCre_ChreCp063';
    'OLD_iCre_ChreCp064'; 'OLD_iCre_ChreCp066'; 'OLD_iCre_ChreCp067';
    'OLD_iCre_ChreCp068'; 'OLD_iCre_ChrepMp01'; 'OLD_iCre_ChrepMp02';
    'OLD_iCre_ChrepMp03'; 'OLD_iCre_ChrepMp04'; 'OLD_iCre_ChrepMp05';
    'OLD_iCre_ChrepMp06'; 'OLD_iCre_ChrepMp07'; 'OLD_iCre_Cre01.g006100.t1.2';
    'OLD_iCre_Cre02.g078300.t1.1'; 'OLD_iCre_Cre02.g095000.t1.2'; 'OLD_iCre_Cre03.g157700.t1.2';
    'OLD_iCre_Cre03.g184400.t1.2'; 'OLD_iCre_Cre06.g278188.t1.1'; 'OLD_iCre_Cre06.g291650.t1.2';
    'OLD_iCre_Cre07.g321050.t1.1'; 'OLD_iCre_Cre07.g333900.t1.1'; 'OLD_iCre_Cre07.g343433.t1.1';
    'OLD_iCre_Cre07.g348200.t1.1'; 'OLD_iCre_Cre09.g390200.t1.1'; 'OLD_iCre_Cre09.g402300.t1.2';
    'OLD_iCre_Cre10.g420350.t1.2'; 'OLD_iCre_Cre10.g420700.t1.2'; 'OLD_iCre_Cre11.g468950.t1.2';
    'OLD_iCre_Cre12.g511200.t1.2'; 'OLD_iCre_Cre12.g546150.t1.2'; 'OLD_iCre_Cre13.g565321.t1.1';
    'OLD_iCre_Cre14.g615550.t1.1'; 'OLD_iCre_Cre16.g650100.t1.2'; 'OLD_iCre_Cre16.g666334.t3.1';
    'OLD_iCre_Cre16.g682050.t1.1'; 'OLD_iCre_Cre16.g682050.t2.1'; 'OLD.iCre.Cre10.g420350.t1.2';
};

model_final = removeGenes(model_final, genesToRemove, false);

printModelStats(model_final, true, true);

[~, textData] = xlsread('ComplementaryData/updated_4414.xlsx','subSystem');

rxnList = textData(2:end,1);
newSubs = textData(2:end,2);

assert(isfield(model_final,'subSystems'), 'Model has no subSystems field');
assert(numel(model_final.subSystems) == numel(model_final.rxns), ...
       'subSystems length does not match rxns');

[tf, loc] = ismember(rxnList, model_final.rxns);

idx_model  = loc(tf);
idx_excel = find(tf);


valid = ~cellfun(@isempty, newSubs(idx_excel));

model_final.subSystems(idx_model(valid)) = newSubs(idx_excel(valid));
for i = 1:numel(model_final.subSystems)

    s = model_final.subSystems{i};

    if ischar(s) || isstring(s)
        model_final.subSystems{i} = {char(s)};

    elseif isempty(s)
        model_final.subSystems{i} = {'Unknown'};

    elseif iscell(s)
        while iscell(s) && numel(s) == 1 && iscell(s{1})
            s = s{1};
        end
        s = s(~cellfun(@isempty, s));

        if isempty(s)
            model_final.subSystems{i} = {'Unknown'};
        else
            model_final.subSystems{i} = cellfun(@char, s, 'UniformOutput', false);
        end
    else
        model_final.subSystems{i} = {'Unknown'};
    end
end

sol = solveLP(model_final,1);
printFluxes(model_final, sol.x, true);
fprintf(['umax = ' num2str(sol.f * -1) ' per hour\n']);


model_final.id= 'iYH2021';
model_final.name = 'Chlamydomonas_reinhardtii_CC4414_GSMM';
exportForGit(model_final,'iYH2021','',{'mat', 'txt', 'xlsx'});
