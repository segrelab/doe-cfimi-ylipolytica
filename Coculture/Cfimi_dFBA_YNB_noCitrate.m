
clear variables
close all
clc

%% Load Data

saveDataName = 'Cfimi_YNB_noCitrate';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% C. fimi Model
load('..\Models\Cfimi_KBase_gfMPIPES_bigg_trsptrs.mat')
model{1} = Cfimi;

% Parameters
load('parameters.mat');

% Medium
load('..\Media\Media_YNB.mat')
medium = arrayfun(@(x) ['EX_' YNB(x).id '_e'],1:numel(YNB),'Uni',false);
media_names = arrayfun(@(x) [YNB(x).id '_e'], 1:length({YNB.id}), 'UniformOutput',false);
media_mM = 1E3.*[YNB(:).M]; % [mM] (1 M * 1E3 mM/1 M)
% C. fimi Essential Metabolites
load('..\Media\essentialMets_Cfimi.mat')
added_names = arrayfun(@(x) [Cfimi_essentialMets(x).id '_e'], 1:length({Cfimi_essentialMets.id}), 'UniformOutput',false);
[~,idx_add] = setdiff(added_names,media_names);
media_names = [media_names added_names(idx_add)];
media_mM = [media_mM, 1E3.*(5E-5.*ones(1,length(added_names(idx_add))))]; % [mM] (1 M * 1E3 mM/1 M)
% Remove Citrate
[~,cit_idx,~] = intersect(media_names,'cit_e');
media_names(cit_idx) = [];
media_mM(cit_idx) = [];
% Cellulose
media_names = [media_names 'cellulose_e'];
media_mM = [media_mM, params.cellulose_mM]; % [mM]

%% Set Exchange Bounds

for numModel = 1:length(model)
    % Find Exchange Reactions
    exch_rxns = identifyExchRxns(model{numModel});
    % Change Bounds
    model{numModel}.lb(exch_rxns) = -100;
    [~,glc_idx] = intersect(model{numModel}.rxns,'EX_glc__D_e');
    model{numModel}.lb(glc_idx) = -10;
    [~,o2_idx] = intersect(model{numModel}.rxns,'EX_o2_e');
    model{numModel}.lb(o2_idx) = -20;
    [~,co2_idx] = intersect(model{1}.rxns,'EX_co2_e');
    model{numModel}.lb(co2_idx) = 0;
end

%% Define Parameter Values

biomass_0 = params.biomass_0.*ones(size(model)); % gCDW
dt = params.dt; % hr
Tend = params.Tend; % hrs
N = Tend/dt;
volume = params.volume; % L
Km = [params.Km.Cfimi]; % mM
Vmax = [params.Vmax.Cfimi]; % mmol/gCDW/hr
max_biomass = params.max_biomass; % gCDW
enzymeModel = 1;

%% Run dFBA

[time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA_cellulase(model,media_names,media_mM.*volume,1,biomass_0,dt,N,volume,Km,Vmax,max_biomass,enzymeModel);
exchMets_mM = exchMets_amt./volume; % [mM] (mmol/L)
exchMets_names = strrep(strrep(erase(exchMets_names,'_e'),'__','-'),'_','-');

save([pwd '\' saveDataName '\' saveDataName '.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');

%% Plots

% Growth
bio_idx = find(model{1}.c);
endGrowth_idx = find(flux{1}(5:end,bio_idx)==0,1)+4;
t = [time(endGrowth_idx),time(end),time(end),time(endGrowth_idx)];

% Exchange Flux & Rename Metabolites to Full Names
exchFlux = cell(size(biomass));
exchFlux_mean = cell(size(biomass));
exchFlux_rxns = cell(size(biomass));
for numModel = 1:length(model)
    % Find Exchange Reactions
    exch_rxns = identifyExchRxns(model{numModel});
    % Exchange Flux
    exchFlux{numModel} = flux{numModel}(:,exch_rxns);
    exchFlux_mean{numModel} = movmean(exchFlux{numModel},3,1); % centered moving average
    exchFlux_mean{numModel}(abs(exchFlux_mean{numModel}) <= 1E-9) = 0;
    exchFlux_rxns{numModel} = model{numModel}.rxnNames(exch_rxns);
    [~,rm_idx,~] = intersect(exchFlux_rxns{numModel},'EX_Biomass_c0');
    exchFlux_rxns{numModel}(rm_idx) = []; exchFlux{numModel}(:,rm_idx) = []; exchFlux_mean{numModel}(:,rm_idx) = [];
    exchFlux_rxns{numModel} = erase(exchFlux_rxns{numModel},' exchange');
    
    % Rename Metabolites to Full Names
    temp = strrep(strrep(erase(erase(model{numModel}.rxns(exch_rxns),'EX_'),'_e'),'__','-'),'_','-');
    [~,rm_idx,~] = intersect(temp,'cpd11416-c0'); % biomass exchange
    temp(rm_idx) = []; exch_rxns(rm_idx) = []; % remove biomass exchange reaction
    [~,mets_idx,temp_idx] = intersect(exchMets_names,temp, 'stable');
    exchMets_names(mets_idx) = erase(model{numModel}.rxnNames(exch_rxns(temp_idx)),' exchange');
end
% Save Additional Data
save([pwd '\' saveDataName '\' saveDataName '.mat'],'exchMets_mM','exchMets_names', 'exchFlux','exchFlux_mean','exchFlux_rxns', '-append')

% Enzyme
[~,cellulase_idx,~] = intersect(exchMets_names,'cellulase');
if ~isempty(cellulase_idx)
    cellulase_mM = exchMets_mM(:,cellulase_idx);
    exchMets_mM(:,cellulase_idx) = [];
    exchMets_names(cellulase_idx) = [];
end

% Metabolite Indices
[~,glc_idx,~] = intersect(exchMets_names,'D-Glucose');
glc_mM = exchMets_mM(:,glc_idx);
[~,cellulose_idx,~] = intersect(exchMets_names,'cellulose');
if ~isempty(cellulose_idx)
    cellulose_mM = exchMets_mM(:,cellulose_idx);
    exchMets_mM(:,cellulose_idx) = [];
    exchMets_names(cellulose_idx) = [];
end

Cfimi_color = [84,39,143]./256;
Ylipolytica_color = [166,54,3]./256;
enzyme_color = [0,109,44]./256;

microbe = {'C. fimi'};
medium = 'BSM+cellulose';
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Biomass Semilog Plot
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass'; ax = axes(fig);
fill(ax,t,[0.01,0.01,30,30],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
semilogy(ax,time,1E3.*biomass{1},'-', 'LineWidth',lineWidth, 'Color',Cfimi_color); hold(ax,'off');
xlim(ax,[0, time(end)])
ax.YScale = 'log'; ylim(ax,[0, 30])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Biomass'],'png')

% Cellulase
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Cellulase'; ax = axes(fig);
fill(ax,t,[0,0,9,9],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
plot(ax,time,cellulase_mM,'-','LineWidth',lineWidth, 'Color',enzyme_color); hold(ax,'off');
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Cellulase Concentration [mM]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Cellulase'],'png')

% Cellulose and Glucose
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Cellulose & Glucose'; ax_c = axes(fig); ax_g = axes(fig);
linkaxes([ax_c,ax_g],'x'); % link x axes
fill(ax_c,t,[0,0,5E-3,5E-3],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax_c,'on');
plot(ax_c,time,cellulose_mM,'k-','LineWidth',lineWidth); hold(ax_c,'off');
plot(ax_g,time,glc_mM,'k--','LineWidth',lineWidth);
xlim(ax_g,[0, time(end)]); ylim(ax_g,[0, 5]);
grid(ax_g,'on'); ax_g.FontSize = axesLabelSize; ax_c.FontSize = axesLabelSize;
ax_g.YAxisLocation = 'right'; ax_g.Color = 'none'; ax_c.XTick = [];
xlabel(ax_g, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax_g, 'Glucose Concentration [mM]', 'FontSize',xyLabelSize)
ylabel(ax_c, 'Cellulose Concentration [mM]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\CelluloseGlucose'],'png')

%% Secreted Metabolites

% Positive Flux
[~,col] = find(exchFlux_mean{1} > 0); col = unique(col);
c = cbrewer('qual','Set1',numel(col));
n = 4;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Positive Flux'; ax = axes(fig); clear h
fill(ax,t,[-10,-10,70,70],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
for ii = 1:numel(col)
    h(ii) = plot(ax,time,exchFlux_mean{1}(:,col(ii)), 'LineWidth',lineWidth, 'Color',c(ii,:));
end; hold(ax,'off');
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Exchange Flux [mmol/(gCDW*hr)]', 'FontSize',xyLabelSize)
lh = legend(h,exchFlux_rxns{1}(col), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\PositiveExchFlux'],'png')

[~,~,met_idx] = intersect(exchFlux_rxns{1}(col),exchMets_names, 'stable');
rm_idx = find(max(exchMets_mM(:,met_idx)) > 100);
rm_metPos = exchMets_names(met_idx(rm_idx));
met_idx(rm_idx) = []; c(rm_idx,:) = [];

% Secreted Metabolites
n = 5;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Secreted Metabolites'; ax = axes(fig); clear h
fill(ax,t,[0,0,1.8,1.8],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
for ii = 1:numel(met_idx)
    h(ii) = plot(ax,time,exchMets_mM(:,met_idx(ii)), 'LineWidth',lineWidth, 'Color',c(ii,:));
end; hold(ax,'off');
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
lh = legend(h,exchMets_names(met_idx), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\SecretedMetabolites'],'png')

%% Consumed Metabolites

% Negative Flux
[~,col] = find(exchFlux_mean{1} < 0); col = unique(col);
c = cbrewer('qual','Set1',numel(col));
n = 6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Negative Flux'; ax = axes(fig); clear h
fill(ax,t,[-10,-10,10,10],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
for ii = 1:numel(col)
    h(ii) = plot(ax,time,exchFlux_mean{1}(:,col(ii)), 'LineWidth',lineWidth, 'Color',c(ii,:));
end; hold(ax,'off');
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Exchange Flux [mmol/(gCDW*hr)]', 'FontSize',xyLabelSize)
lh = legend(h,exchFlux_rxns{1}(col), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\NegativeExchFlux'],'png')

[~,~,met_idx] = intersect(exchFlux_rxns{1}(col),exchMets_names, 'stable');
rm_idx = find(max(exchMets_mM(:,met_idx)) > 100);
rm_metNeg = exchMets_names(met_idx(rm_idx));
met_idx(rm_idx) = [];

% Consumed Metabolites - 1
idx = find(max(exchMets_mM(:,met_idx)) > 0.1);
n = 7;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Consumed Metabolites'; ax = axes(fig); clear h
fill(ax,t,[0,0,60,60],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax,'on');
for ii = 1:numel(idx)
    h(ii) = plot(ax,time,exchMets_mM(:,met_idx(idx(ii))), 'LineWidth',lineWidth, 'Color',c(idx(ii),:));
end; hold(ax,'off');
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
lh = legend(h,exchMets_names(met_idx(idx)), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\ConsumedMetabolites1'],'png')

% Consumed Metabolites - 2
idx_mM = find(max(exchMets_mM(:,met_idx)) <= 0.1 & max(exchMets_mM(:,met_idx)) >= 1E-3);
idx_uM = find(max(exchMets_mM(:,met_idx)) < 1E-3);
n = 8;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Consumed Metabolites'; ax_mM = axes(fig); ax_uM = axes(fig); clear h
linkaxes([ax_mM,ax_uM],'x'); % link x axes
fill(ax_mM,t,[0,0,1,1],0.5.*[1,1,1], 'LineStyle','none', 'FaceAlpha',0.25); hold(ax_mM,'on');
for ii = 1:numel(idx_mM)
    h(ii) = plot(ax_mM,time,exchMets_mM(:,met_idx(idx_mM(ii))), 'LineWidth',lineWidth, 'Color',c(idx_mM(ii),:));
end; hold(ax_mM,'off');
for ii = 1:numel(idx_uM)
    g(ii) = plot(ax_uM,time,1E3.*exchMets_mM(:,met_idx(idx_uM(ii))), 'LineWidth',lineWidth, 'Color',c(idx_uM(ii),:)); hold(ax_uM,'on');
end; hold(ax_uM,'off');
xlim(ax_uM,[0, time(end)]); ylim(ax_mM,[0,0.05]); ylim(ax_uM,[0,1])
grid(ax_uM,'on'); ax_uM.FontSize = axesLabelSize; ax_mM.FontSize = axesLabelSize;
ax_uM.YAxisLocation = 'right'; ax_uM.Color = 'none'; ax_mM.XTick = [];
xlabel(ax_uM, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax_mM, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
ylabel(ax_uM, 'Metabolite Concentration [\muM]', 'FontSize',xyLabelSize)
lh = legend([h,g],exchMets_names(met_idx([idx_mM, idx_uM])), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\ConsumedMetabolites2'],'png')




