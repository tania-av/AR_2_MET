%% FASTCORMICS RNA-seq
% _(c) Dr. Maria Pires Pacheco 2016_
% 
% _Example script adapted by Tamara Bintener_
%% Disclaimer
% In order to run FASTCORMICS RNA-seq you will need
%% 
% * Matlab
% * Compatible Cplex version, added to the Matlab path
% * Curve fitting toolbox
% * Cobra Toolbox installed (https://opencobra.github.io/cobratoolbox/latest/installation.html)
%% 
% The example data was downloaded from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009 
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009> (FPKM values 
% of normal tissue). The first two samples (TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07) 
% will serve as an example was already imported, split into 3 variables, and saved 
% under _example_FPKM.mat ._
%% 
% * colnames:   cell array with the sample names (here TCGA06067511A32RA36H07)
% * rownames:  cell aray with the gene IDs
% * fpkm:           fpkm values for the samples, size(fpkm) = lenghth(geneIDs) 
% length(colnames)

%% Setup
clear all
close all
initCobraToolbox
addpath(genpath(pwd)) % add all subfolders to path

%% load the data

% %TRICLOSAN
% b = readtable('data1.csv');
% Genes
% rownames = b.Names;
% colnames = {'NAB' 'CAB_0.02' 'CAB_2.0'} %0, 0.02, 2.00 mg/L
% fpkm = [b.NAB b.CAB_1 b.CAB_2]

%C. CAMPHORA
b = readtable('data2.csv');
%Genes
rownames = b.Names;
colnames = {'NAB' 'CAB_1' 'CAB_2' 'CAB_3' 'CAB_4'} %0, 0.02, 2.00 mg/L
fpkm = [b.NAB b.CAB_1 b.CAB_2 b.CAB_3 b.CAB_4]

% %TAT
% b = readtable('data4.csv');
% %Genes
% rownames = b.Names;
% colnames = {'NAB_M9' 'CAB_M9' 'NAB_M9T' 'CAB_M9T'} %0, 1.56mg/mL, 0, 1.56mg/mL
% fpkm = [b.NAB_M9 b.CAB_M9 b.NAB_M9T b.CAB_M9T]

figflag = 1; % set figure flag: 1 to output and save density plot
%% Data discretization

discretized = discretize_FPKM(fpkm, colnames) % no figures

%% FASTCORMICS
%% Prepare FASTCORMICS
% needed:
%% 
% Model:
Cmodel_original = readCbModel('iML1515.xml');
Cmodel_original.rev = Cmodel_original.lb<1;
sol_original = optimizeCbModel(Cmodel_original)

Cmodel_original = changeObjective(Cmodel_original, 'BIOMASS_Ec_iML1515_core_75p37M')
optimizeCbModel(Cmodel_original)

%From rFASTCORMICS_V-2020 example
A = fastcc_4_rfastcormics(Cmodel_original, 1e-10,0) % create consistent model by running FASTCC (Vlassis et al., 2014)
Cmodel = removeRxns(Cmodel_original, Cmodel_original.rxns(setdiff(1:numel(Cmodel_original.rxns),A)))
optimizeCbModel(Cmodel)

%%
dico = cell2table([rownames rownames]); % dictionary to map the rownname identifier to the genes in the model
epsilon = 1e-5;
already_mapped_tag = 0;
consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples
%%  Set optional settings such as:

optional_settings.func = {'BIOMASS_Ec_iML1515_core_75p37M'}; % forced reactions
biomass_rxn = 'BIOMASS_Ec_iML1515_core_75p37M';

%% Create models

% Multiply the bounds by 10 to minimize numerical issues:
Cmodel_10 = Cmodel;
Cmodel_10.lb = Cmodel_10.lb*10;
Cmodel_10.ub = Cmodel_10.ub*10;

% single models:
for i = 1:numel(colnames) %for each sample
    [model_out{i}, A_keep{i}] = fastcormics_RNAseq(Cmodel_10, discretized(:,i), ...
        rownames, dico , biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
end

%% Results

% % TRICLOSAN
% length(model_out{1}.rxns)
% length(model_out{2}.rxns)
% length(model_out{3}.rxns)
% 
% nab = model_out{1};
% nab.lb = nab.lb/10;
% nab.ub = nab.ub/10;
% optimizeCbModel(nab)
% 
% cab1 = model_out{2};
% cab1.lb = cab1.lb/10;
% cab1.ub = cab1.ub/10;
% optimizeCbModel(cab1)
% 
% cab2 = model_out{3};
% cab2.lb = cab2.lb/10;
% cab2.ub = cab2.ub/10;
% optimizeCbModel(cab2)

% %C. CAMPHORA
% length(model_out{1}.rxns)
% length(model_out{2}.rxns)
% length(model_out{3}.rxns)
% length(model_out{4}.rxns)
% length(model_out{5}.rxns)
% 
% nab = model_out{1};
% nab.lb = nab.lb/10;
% nab.ub = nab.ub/10;
% optimizeCbModel(nab)
% 
% cab1 = model_out{2};
% cab1.lb = cab1.lb/10;
% cab1.ub = cab1.ub/10;
% optimizeCbModel(cab1)
% 
% cab2 = model_out{3};
% cab2.lb = cab2.lb/10;
% cab2.ub = cab2.ub/10;
% optimizeCbModel(cab2)
% 
% cab3 = model_out{4};
% cab3.lb = cab3.lb/10;
% cab3.ub = cab3.ub/10;
% optimizeCbModel(cab3)
% 
% cab4 = model_out{5};
% cab4.lb = cab4.lb/10;
% cab4.ub = cab4.ub/10;
% optimizeCbModel(cab4)

% TAT
length(model_out{1}.rxns)
length(model_out{2}.rxns)
length(model_out{3}.rxns)
length(model_out{4}.rxns)

length(model_out{1}.genes)
length(model_out{2}.genes)
length(model_out{3}.genes)
length(model_out{4}.genes)

nab = model_out{1};
nab.lb = nab.lb/10;
nab.ub = nab.ub/10;
optimizeCbModel(nab)

cab1 = model_out{2};
cab1.lb = cab1.lb/10;
cab1.ub = cab1.ub/10;
optimizeCbModel(cab1)

nab2 = model_out{3};
nab2.lb = nab2.lb/10;
nab2.ub = nab2.ub/10;
optimizeCbModel(nab2)

cab2 = model_out{4};
cab2.lb = cab2.lb/10;
cab2.ub = cab2.ub/10;
optimizeCbModel(cab2)

% generic models:

[model_out_consensus, A_keep_consensus] = fastcormics_RNAseq(Cmodel, discretized, ...
    rownames, dico , biomass_rxn, already_mapped_tag, consensus_proportion, epsilon);
%%

delete clone*.log %delete some files generated by cplex

