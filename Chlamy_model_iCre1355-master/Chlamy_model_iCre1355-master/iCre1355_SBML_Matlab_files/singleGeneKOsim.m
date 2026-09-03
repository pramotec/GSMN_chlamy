%
% This script performs a single gene analysis of the given models via 
% singleGeneDeletion method of the COBRA toolbox
%
% PREREQUISITES:
%   - have COBRA toolbox installed and initialized (run initCobraToolbox)
%   - make sure the files iCre1355_{auto|mixo|hetero}.xml as well as 
%     transcript_gene_map.mat are in your current folder
%
% Output description:
%
% each singleKOsim.{auto|mixo|hetero} contains the following data (see 
% also singleGeneDeletion.m from the COBRA toolbox)
%
% grRatio       Computed growth rate ratio between deletion strain and wild type
% grRateKO      Deletion strain growth rates (1/h)
% grRateWT      Wild type growth rate (1/h)
% hasEffect     Does a gene deletion affect anything (i.e. are any reactions
%               removed from the model)
% delRxns       List of deleted reactions for each gene KO
% fluxSolution  FBA/MOMA/lMOMA fluxes for KO strains

%
model_auto   = readCbModel('iCre1355_auto.xml');
model_mixo   = readCbModel('iCre1355_mixo.xml');
model_hetero = readCbModel('iCre1355_hetero.xml');
load transcript_gene_map

model_auto.c(2392) = 1;
model_mixo.c(2393) = 1;
model_hetero.c(2394) = 1;

% process list of KO candidate genes via transcript_gene_map
tmp = unique(transcript_gene_map(2:end,2));
geneKOlist = cell(length(tmp),1);
for i=1:length(tmp)
    geneKOlist{i} = transcript_gene_map(...
        strcmp(tmp(i,1),transcript_gene_map(:,2)),1);
end
clear tmp

[singleKOsim.auto.grRatio,...
    singleKOsim.auto.grRateKO,...
    singleKOsim.auto.grRateWT,...
    singleKOsim.auto.hasEffect,...
    singleKOsim.auto.delRxns,...
    singleKOsim.auto.fluxSolution] = singleGeneDeletion(model_auto, 'FBA', geneKOlist);
[singleKOsim.mixo.grRatio,...
    singleKOsim.mixo.grRateKO,...
    singleKOsim.mixo.grRateWT,...
    singleKOsim.mixo.hasEffect,...
    singleKOsim.mixo.delRxns,...
    singleKOsim.mixo.fluxSolution] = singleGeneDeletion(model_mixo, 'FBA', geneKOlist);
[singleKOsim.hetero.grRatio,...
    singleKOsim.hetero.grRateKO,...
    singleKOsim.hetero.grRateWT,...
    singleKOsim.hetero.hasEffect,...
    singleKOsim.hetero.delRxns,...
    singleKOsim.hetero.fluxSolution] = singleGeneDeletion(model_hetero, 'FBA', geneKOlist);