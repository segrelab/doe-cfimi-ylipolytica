
clear variables
close all
clc

%% Load Data

% FBA Simulations
bsm_fba = load('Cfimi_BiologFBA_BSM.mat');
mpipes_fba = load('Cfimi_BiologFBA_MPIPES.mat');

% Experiments
bsm_exp = load('Cfimi_BiologExp_BSM.mat');
mpipes_exp = load('Cfimi_BiologExp_MPIPES.mat');

% Calculations
expPos = @(media) [sum(arrayfun(@(x) ismember(media{x,2},'+', 'rows') && ismember(media{x,3},'+', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'-', 'rows') && ismember(media{x,3},'+', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'n/a', 'rows') && ismember(media{x,3},'+', 'rows'),1:size(media,1)))];
expNeg = @(media) [sum(arrayfun(@(x) ismember(media{x,2},'+', 'rows') && ismember(media{x,3},'-', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'-', 'rows') && ismember(media{x,3},'-', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'n/a', 'rows') && ismember(media{x,3},'-', 'rows'),1:size(media,1)))];
expN_A = @(media) [sum(arrayfun(@(x) ismember(media{x,2},'+', 'rows') && ismember(media{x,3},'n/a', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'-', 'rows') && ismember(media{x,3},'n/a', 'rows'),1:size(media,1)));
    sum(arrayfun(@(x) ismember(media{x,2},'n/a', 'rows') && ismember(media{x,3},'n/a', 'rows'),1:size(media,1)))];

%% Confusion Matrices

% BSM
[~,idx_fba,idx_exp] = intersect(bsm_fba.carbon_name,bsm_exp.carbon_name);
bsm = [bsm_fba.carbon_name(idx_fba)', bsm_fba.growth_class(idx_fba), bsm_exp.growth_class(idx_exp)];
T_bsm = table(expPos(bsm),expNeg(bsm),expN_A(bsm));
T_bsm.Properties.VariableNames = {'expPos','expNeg','expN_A'};
T_bsm.Properties.RowNames = {'fbaPos','fbaNeg','fbaN_A'};
T_bsm

% MPIPES
[~,idx_fba,idx_exp] = intersect(mpipes_fba.carbon_name,mpipes_exp.carbon_name);
mpipes = [mpipes_fba.carbon_name(idx_fba)', mpipes_fba.growth_class(idx_fba), mpipes_exp.growth_class(idx_exp)];
T_mpipes = table(expPos(mpipes),expNeg(mpipes),expN_A(mpipes));
T_mpipes.Properties.VariableNames = {'expPos','expNeg','expN_A'};
T_mpipes.Properties.RowNames = {'fbaPos','fbaNeg','fbaN_A'};
T_mpipes

