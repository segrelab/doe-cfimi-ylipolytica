
clear variables
close all
clc

%% Load Data

% Model
load('..\Models\Cfimi_KBase_gfMPIPES_bigg_trsptrs.mat')

% Essential Metabolites
load('..\Media\essentialMets_Cfimi.mat')

% Biolog Sources
load('..\Media\Biolog_CNPS_sources.mat')
carbon_id = arrayfun(@(x) ['EX_' CNPS_sources.carbon(x).bigg_id '_e'],1:96,'Uni',false);
carbon_name = arrayfun(@(x) CNPS_sources.carbon(x).chemical,1:96,'Uni',false);
carbon_formula = arrayfun(@(x) CNPS_sources.carbon(x).formula,1:96,'Uni',false);

% Medium
load('..\Media\Media_modifiedBSM.mat')
medium = arrayfun(@(x) ['EX_' modified_BSM(x).id '_e'],1:numel(modified_BSM),'Uni',false);

%% Set Default Exchange Reaction Lower Bounds

numC = zeros(size(carbon_formula));
for ii = 2:96
    sol = strsplit(carbon_formula{ii},'C');
    sol = strsplit(sol{2},'H');
    if isempty(str2num(sol{1}))
        numC(ii) = 1;
    else
        numC(ii) = str2num(sol{1});
    end
end
clear temp

% Exchange Reactions
[exch_rxns,~] = identifyExchRxns(Cfimi);
Cfimi.lb(exch_rxns) = 0;

% Find Reactions for Essential Metabolites
essential_rxns = arrayfun(@(x) ['EX_' Cfimi_essentialMets(x).id '_e'],1:numel(Cfimi_essentialMets),'Uni',false);
[~,idx,~] = intersect(Cfimi.rxns,essential_rxns);
Cfimi.lb(idx) = -1000;

% Medium
[~,idx,~] = intersect(Cfimi.rxns,medium);
Cfimi.lb(idx) = -1000;
[~,idx_O2,~] = intersect(Cfimi.rxns,'EX_o2_e');
Cfimi.lb(idx_O2) = -20;

clear essential_rxns idx*

%% Carbon Sources

% Pre-Allocate
growth_rate = NaN(numel(carbon_id),1);
growth_class = cell(numel(carbon_id),1);

% Negative Control
model = Cfimi;
sol = FBA(model,[],true);
growth_rate(1) = sol.objectiveValue;
if sol.objectiveValue==0 || isnan(sol.objectiveValue)
    growth_class{1} = '-';
elseif sol.objectiveValue>0
    growth_class{1} = '+';
end
% Carbon Sources
for ii = 2:numel(carbon_id)
    model = Cfimi;
    [~,idx_C,~] = intersect(Cfimi.rxns,carbon_id{ii}); % find carbon source
    if isempty(idx_C)
        growth_rate(ii) = 0;
        growth_class{ii} = 'n/a';
    else
        model.lb(idx_C) = -10.*(numC(ii)/6); % set lower bound

        sol = FBA(model,[],true);
        growth_rate(ii) = sol.objectiveValue;
        if sol.objectiveValue==0 || isnan(sol.objectiveValue)
            growth_class{ii} = '-';
        elseif sol.objectiveValue>0
            growth_class{ii} = '+';
        end
    end
end

save([pwd '\Cfimi_BiologFBA_BSM.mat'],'carbon_name','growth_rate','growth_class')






