
clear variables
close all
clc

%% Load Data

saveDataName = 'Cfimi_Ylipolytica_YNB';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% C. fimi Model
load('..\Models\Cfimi_KBase_gfMPIPES_bigg_trsptrs.mat')
model{1} = Cfimi;
% Y. lipolytica Model
load('..\Models\Ylipolytica_iMK735_bigg.mat')
model{2} = Ylipolytica;

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
% Y. lipolytica Essential Metabolites
load('..\Media\essentialMets_Ylipolytica.mat')
added_names = arrayfun(@(x) [Ylipolytica_essentialMets(x).id '_e'], 1:length({Ylipolytica_essentialMets.id}), 'UniformOutput',false);
[~,idx_add] = setdiff(added_names,media_names);
media_names = [media_names added_names(idx_add)];
media_mM = [media_mM, 1E3.*(5E-5.*ones(1,length(added_names(idx_add))))]; % [mM] (1 M * 1E3 mM/1 M)
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
    [~,co2_idx] = intersect(model{numModel}.rxns,'EX_co2_e');
    model{numModel}.lb(co2_idx) = 0;
end

%% Define Parameter Values

biomass_0 = params.biomass_0.*ones(size(model)); % gCDW
dt = params.dt; % hr
Tend = params.Tend; % hrs
N = Tend/dt;
volume = params.volume; % L
Km = [params.Km.Cfimi, params.Km.Ylipolytica]; % mM
Vmax = [params.Vmax.Cfimi, params.Vmax.Ylipolytica]; % mmol/gCDW/hr
max_biomass = params.max_biomass; % gCDW
enzymeModel = 1;

%% Run dFBA

[time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA_cellulase(model,media_names,media_mM.*volume,1,biomass_0,dt,N,volume,Km,Vmax,max_biomass,enzymeModel);
exchMets_mM = exchMets_amt./volume; % [mM] (mmol/L)
exchMets_names = strrep(strrep(erase(exchMets_names,'_e'),'__','-'),'_','-');

save([pwd '\' saveDataName '\' saveDataName '.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');

%% Plots

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

microbe = {'C. fimi'; 'Y. lipolytica'};
medium = 'YNB+cellulose';
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Biomass Semilog Plot
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass'; ax = axes(fig);
semilogy(ax,time,1E3.*biomass{1},'-', 'LineWidth',lineWidth, 'Color',Cfimi_color); hold(ax, 'on');
semilogy(ax,time,1E3.*biomass{2},'--', 'LineWidth',lineWidth, 'Color',Ylipolytica_color); hold(ax, 'off');
xlim(ax,[0, time(end)])
ylim(ax,[0, 3E3])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
lh = legend(microbe, 'Location','Best'); lh.Box = 'Off';
title(ax, medium, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\Biomass'],'svg')

% Cellulase
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Cellulase'; ax = axes(fig);
plot(ax,time,cellulase_mM,'-','LineWidth',lineWidth, 'Color',enzyme_color)
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Cellulase Concentration [mM]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Cellulase'],'svg')

% Cellulose and Glucose
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Cellulose & Glucose'; ax_c = axes(fig); ax_g = axes(fig);
linkaxes([ax_c,ax_g],'x'); % link x axes
plot(ax_c,time,cellulose_mM,'k-','LineWidth',lineWidth);
plot(ax_g,time,glc_mM,'k--','LineWidth',lineWidth);
xlim(ax_g,[0, time(end)]); ylim(ax_g,[0,1])
grid(ax_g,'on'); ax_g.FontSize = axesLabelSize; ax_c.FontSize = axesLabelSize;
ax_g.YAxisLocation = 'right'; ax_g.Color = 'none'; ax_c.XTick = [];
xlabel(ax_g, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax_g, 'Glucose Concentration [mM]', 'FontSize',xyLabelSize)
ylabel(ax_c, 'Cellulose Concentration [mM]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\CelluloseGlucose'],'svg')

%% Secreted Metabolites

% Positive Flux
n = 4; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Positive Flux';
for numModel = 1:numel(model)
    [~,col{numModel}] = find(exchFlux_mean{numModel} > 0); col{numModel} = unique(col{numModel});
    c = cbrewer('qual','Set1',numel(col{numModel}));
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(col{numModel})
        h(ii) = plot(ax,time,exchFlux_mean{numModel}(:,col{numModel}(ii)), 'LineWidth',f*lineWidth, 'Color',c(ii,:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Exchange Flux [mmol/(gCDW*hr)]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchFlux_rxns{numModel}(col{numModel}), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\PositiveExchFlux'],'svg')

% Secreted Metabolites
n = 5; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Secreted Metabolites';
for numModel = 1:numel(model)
    c = cbrewer('qual','Set1',numel(col{numModel}));
    [~,~,met_idx{numModel}] = intersect(exchFlux_rxns{numModel}(col{numModel}),exchMets_names, 'stable');
    rm_idx = find(max(exchMets_mM(:,met_idx{numModel})) > 100);
    rm_metPos{numModel} = exchMets_names(met_idx{numModel}(rm_idx));
    met_idx{numModel}(rm_idx) = []; c(rm_idx,:) = [];
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(met_idx{numModel})
        h(ii) = plot(ax,time,exchMets_mM(:,met_idx{numModel}(ii)), 'LineWidth',f*lineWidth, 'Color',c(ii,:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchMets_names(met_idx{numModel}), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\SecretedMetabolites'],'svg')

%% Consumed Metabolites

% Negative Flux
n = 6; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Negative Flux'; clear c
for numModel = 1:numel(model)
    [~,col{numModel}] = find(exchFlux_mean{numModel} < 0); col{numModel} = unique(col{numModel});
    c{numModel} = cbrewer('qual','Set1',numel(col{numModel}));
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(col{numModel})
        h(ii) = plot(ax,time,exchFlux_mean{numModel}(:,col{numModel}(ii)), 'LineWidth',f*lineWidth, 'Color',c{numModel}(ii,:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Exchange Flux [mmol/(gCDW*hr)]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchFlux_rxns{numModel}(col{numModel}), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\NegativeExchFlux'],'svg')

% Consumed Metabolites - 1
n = 7; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Consumed Metabolites';
for numModel = 1:numel(model)
    [~,~,met_idx{numModel}] = intersect(exchFlux_rxns{numModel}(col{numModel}),exchMets_names, 'stable');
    rm_idx = find(max(exchMets_mM(:,met_idx{numModel})) > 100);
    rm_metNeg{numModel} = exchMets_names(met_idx{numModel}(rm_idx));
    met_idx{numModel}(rm_idx) = []; c{numModel}(rm_idx,:) = [];
    
    idx{numModel} = find(max(exchMets_mM(:,met_idx{numModel})) > 0.1);
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(idx{numModel})
        h(ii) = plot(ax,time,exchMets_mM(:,met_idx{numModel}(idx{numModel}(ii))), 'LineWidth',f*lineWidth, 'Color',c{numModel}(idx{numModel}(ii),:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchMets_names(met_idx{numModel}(idx{numModel})), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\ConsumedMetabolites1'],'svg')

% Consumed Metabolites - 2
n = 8; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Consumed Metabolites';
for numModel = 1:numel(model)
    idx{numModel} = find(max(exchMets_mM(:,met_idx{numModel})) <= 0.1 & max(exchMets_mM(:,met_idx{numModel})) >= 1E-3);
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(idx{numModel})
        h(ii) = plot(ax,time,exchMets_mM(:,met_idx{numModel}(idx{numModel}(ii))), 'LineWidth',f*lineWidth, 'Color',c{numModel}(idx{numModel}(ii),:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchMets_names(met_idx{numModel}(idx{numModel})), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\ConsumedMetabolites2'],'svg')

% Consumed Metabolites - 3
n = 9; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Consumed Metabolites';
for numModel = 1:numel(model)
    idx{numModel} = find(max(exchMets_mM(:,met_idx{numModel})) < 1E-3);
    ax = subplot(numel(model),1,numModel); clear h
    for ii = 1:numel(idx{numModel})
        h(ii) = plot(ax,time,exchMets_mM(:,met_idx{numModel}(idx{numModel}(ii))), 'LineWidth',f*lineWidth, 'Color',c{numModel}(idx{numModel}(ii),:)); hold(ax,'on');
    end; hold(ax,'off');
    xlim(ax,[0, time(end)])
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
    ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',f*xyLabelSize)
    title(ax, microbe{numModel}, 'FontSize',f*titleSize)
    lh = legend(h,exchMets_names(met_idx{numModel}(idx{numModel})), 'Location','Best'); lh.Box = 'Off';
end
saveas(fig,[pwd '\' saveDataName '\ConsumedMetabolites3'],'svg')

%% Exchanged Metabolites

markerSize = 15;
arrowHead = 4;

[rxns,idx1,idx2] = intersect(exchFlux_rxns{1},exchFlux_rxns{2});
rm_idx = find(sum(abs(exchFlux_mean{1}(:,idx1)))==0 & sum(abs(exchFlux_mean{2}(:,idx2)))==0);
rxns(rm_idx) = []; idx1(rm_idx) = []; idx2(rm_idx) = [];

fluxType = zeros(N+1,numel(rxns));

% Crossfeedogram
n = 10; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Crossfeedogram';
for r = 1:numel(rxns)
    ax = subplot(5,3,r);
    [~,~,l,~,fluxType(:,r)] = crossfeedogram(exchFlux_mean{1}(:,idx1(r)),exchFlux_mean{2}(:,idx2(r)),5,fig,ax,markerSize,arrowHead,f*lineWidth);
    l.LineWidth = 0.5;
    grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    ax.XLabel.String = []; ax.YLabel.String = [];
    title(ax,rxns(r), 'FontSize',f*titleSize)
end
[~,hx] = suplabel([microbe{1} ' Exchange Flux']); set(hx, 'FontSize',xyLabelSize);
[~,hy] = suplabel([microbe{2} ' Exchange Flux'],'y'); set(hy, 'FontSize',xyLabelSize);
saveas(fig,[pwd '\' saveDataName '\Crossfeedogram'],'svg')

if isdir([pwd '\' saveDataName '\Metabolites'])==0; mkdir([pwd '\' saveDataName '\Metabolites']); end
cmap = [[233,233,233];[0,69,41];[21,122,62];[81,181,102];[8,48,107];[22,102,174];[79,156,204];[254,190,72];[208,112,12]]./256;
% Metabolite Concentration
for r = 1:numel(rxns)
    [~,~,idx] = intersect(rxns(r),exchMets_names);
    n = 10+r;
    if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
    fig = figure(n); fig.Name = rxns{r};
    ax1 = axes(fig); ax2 = axes(fig);
    linkaxes([ax1,ax2],'x'); % link x axes
    imagesc(ax1,time,ones(N+1,1),fluxType(:,r)')
    colormap(ax1,cmap); caxis(ax1,[0, 8]); cbar = colorbar(ax1);
    cbar.Label.String = 'Exchange Type'; cbar.Ticks = 0.45:8/9:8;
    cbar.TickLabels = {'\rightarrow','Cf,Yl\rightarrow','Cf\rightarrow','Yl\rightarrow',...
        '\rightarrowCf,Yl','\rightarrowCf','\rightarrowYl','Yl\rightarrowCf','Cf\rightarrowYl'};
    plot(ax2,time,exchMets_mM(:,idx),'k', 'LineWidth',lineWidth)
    xlim(ax1,[0, time(end)])
    ax2.Position = ax1.Position;
    grid(ax1,'on'); ax1.FontSize = axesLabelSize; ax2.FontSize = axesLabelSize;
    ax1.Color = 'none'; ax2.Color = 'none'; ax1.YTick = []; ax2.XTick = [];
    xlabel(ax1, 'Time [hours]', 'FontSize',xyLabelSize)
    ylabel(ax2, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
    title(ax1,rxns(r), 'FontSize',titleSize)
    saveas(fig,[pwd '\' saveDataName '\Metabolites\' rxns{r}],'svg')
end


