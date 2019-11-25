
clear variables
close all
clc

%% Load Data

% FBA Simulations
mpipes_fba = load('Ylipolytica_BiologFBA_MPIPES.mat');
mpipes_noCitrate_fba = load('Ylipolytica_BiologFBA_MPIPES_noCitrate.mat');

% Experiments
mpipes_exp = load('Ylipolytica_BiologExp_MPIPES.mat');

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

% MPIPES
[~,idx_fba,idx_exp] = intersect(mpipes_fba.carbon_name,mpipes_exp.carbon_name);
mpipes = [mpipes_fba.carbon_name(idx_fba)', mpipes_fba.growth_class(idx_fba), mpipes_exp.growth_class(idx_exp)];
T_mpipes = table(expPos(mpipes),expNeg(mpipes),expN_A(mpipes));
T_mpipes.Properties.VariableNames = {'expPos','expNeg','expN_A'};
T_mpipes.Properties.RowNames = {'fbaPos','fbaNeg','fbaN_A'};
T_mpipes

% MPIPES
[~,idx_fba,idx_exp] = intersect(mpipes_noCitrate_fba.carbon_name,mpipes_exp.carbon_name);
mpipes_noCit = [mpipes_noCitrate_fba.carbon_name(idx_fba)', mpipes_noCitrate_fba.growth_class(idx_fba), mpipes_exp.growth_class(idx_exp)];
T_mpipes_noCit = table(expPos(mpipes_noCit),expNeg(mpipes_noCit),expN_A(mpipes_noCit));
T_mpipes_noCit.Properties.VariableNames = {'expPos','expNeg','expN_A'};
T_mpipes_noCit.Properties.RowNames = {'fbaPos','fbaNeg','fbaN_A'};
T_mpipes_noCit

