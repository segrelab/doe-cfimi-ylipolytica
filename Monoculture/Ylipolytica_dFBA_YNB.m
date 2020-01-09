
clear variables
close all
clc

%% Load Data

saveDataName = 'Ylipolytica_YNB';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% Model
load('..\Models\Ylipolytica_iMK735_bigg.mat')
model{1} = Ylipolytica;

% Parameters
load('parameters.mat');

% Medium
load('..\Media\Media_YNB.mat')
medium = arrayfun(@(x) ['EX_' YNB(x).id '_e'],1:numel(YNB),'Uni',false);
media_names = arrayfun(@(x) [YNB(x).id '_e'], 1:length({YNB.id}), 'UniformOutput',false);
media_mM = 1E3.*[YNB(:).M]; % [mM] (1 M * 1E3 mM/1 M)
% Essential Metabolites
load('..\Media\essentialMets_Ylipolytica.mat')
added_names = arrayfun(@(x) [Ylipolytica_essentialMets(x).id '_e'], 1:length({Ylipolytica_essentialMets.id}), 'UniformOutput',false);
[~,idx_add] = setdiff(added_names,media_names);
media_names = [media_names added_names(idx_add)];
media_mM = [media_mM, 1E3.*(5E-5.*ones(1,length(added_names(idx_add))))]; % [mM] (1 M * 1E3 mM/1 M)
% Sugars
% sugar_names = {'glc__D';'cellb';{'arab__L';'arab__D'};'man';'gal';'xyl__D'}; % glucose, cellobiose, arabinose, mannose, galactose, xylose
sugar_names = {'glc__D';{'arab__L';'arab__D'};'man';'gal';'xyl__D'}; % glucose, arabinose, mannose, galactose, xylose
sugar_mM = 5.5; % [mM];

%% Set Exchange Bounds

for numModel = 1:length(model)
    % Find Exchange Reactions
    exch_rxns = identifyExchRxns(model{numModel});
    % Change Bounds
    model{numModel}.lb(exch_rxns) = -100;
    [~,o2_idx] = intersect(model{numModel}.rxns,'EX_o2_e');
    model{numModel}.lb(o2_idx) = -20;
    [~,co2_idx] = intersect(model{numModel}.rxns,'EX_co2_e');
    model{numModel}.lb(co2_idx) = 0;
end

%% Define Parameter Values

biomass_0 = params.biomass_0.*ones(size(model)); % gCDW
dt = params.dt; % hr
Tend = params.Tend; % hrs
N = Tend/dt;
volume = params.volume; % L
Km = params.Km.Ylipolytica; % mM
Vmax = params.Vmax.Ylipolytica; % mmol/gCDW/hr
max_biomass = params.max_biomass; % gCDW
enzymeModel = 0;

%% Run dFBA - Individual Sugars

base_media = media_names;
for s = 1:numel(sugar_names)
    media_names = base_media;
    if size(sugar_names{s},1)==1
        [~,c_idx] = intersect(model{numModel}.rxns,['EX_' sugar_names{s} '_e']);
        media_names = [media_names [sugar_names{s} '_e']];
        media_mM = [media_mM, sugar_mM];
    else
        for ii = 1:numel(sugar_names{s})
            [~,c_idx(ii)] = intersect(model{numModel}.rxns,['EX_' sugar_names{s}{ii} '_e']);
            media_names = [media_names [sugar_names{s}{ii} '_e']];
            media_mM = [media_mM, sugar_mM];
        end
    end
    model{numModel}.lb(c_idx) = -10;

    [time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA_cellulase(model,media_names,media_mM.*volume,1,biomass_0,dt,N,volume,Km,Vmax,max_biomass,enzymeModel);
    exchMets_mM = exchMets_amt./volume; % [mM] (mmol/L)
    exchMets_names = strrep(strrep(erase(exchMets_names,'_e'),'__','-'),'_','-');
    
    if size(sugar_names{s},1)==1
        save([pwd '\' saveDataName '\' saveDataName '_' sugar_names{s} '.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');
    else
        temp = strsplit(sugar_names{s}{1},'_');
        save([pwd '\' saveDataName '\' saveDataName '_' temp{1} '.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');
    end
    model{numModel}.lb(c_idx) = -100; clear c_idx
end

%% Run dFBA - All Sugars

media_names = base_media;
c_idx = [];
for s = 1:numel(sugar_names)
    if size(sugar_names{s},1)==1
        [~,temp] = intersect(model{numModel}.rxns,['EX_' sugar_names{s} '_e']);
        c_idx = [c_idx; temp];
        media_names = [media_names [sugar_names{s} '_e']];
        media_mM = [media_mM, sugar_mM];
    else
        for ii = 1:numel(sugar_names{s})
            [~,temp] = intersect(model{numModel}.rxns,['EX_' sugar_names{s}{ii} '_e']);
            c_idx = [c_idx; temp];
            media_names = [media_names [sugar_names{s}{ii} '_e']];
            media_mM = [media_mM, sugar_mM];
        end
    end
    model{numModel}.lb(c_idx) = -10;

    [time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA_cellulase(model,media_names,media_mM.*volume,1,biomass_0,dt,N,volume,Km,Vmax,max_biomass,enzymeModel);
    exchMets_mM = exchMets_amt./volume; % [mM] (mmol/L)
    exchMets_names = strrep(strrep(erase(exchMets_names,'_e'),'__','-'),'_','-');
    
    save([pwd '\' saveDataName '\' saveDataName '_all.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');
end









