
clear variables
close all
clc

%% Load Data

saveDataName = 'Cfimi';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% Experimental Data
exp_bsm = load([saveDataName '_BSM_exp.mat']);
exp_ynb = load([saveDataName '_YNB_exp.mat']);

% Simulation Data
data = load([saveDataName '_YNB/' saveDataName '_YNB.mat']);
time = data.time; biomass{1} = data.biomass{1};
data = load([saveDataName '_BSM/' saveDataName '_BSM.mat']);
biomass{2} = data.biomass{1};

%% Plots

microbe = {'C. fimi'};
medium = {'YNB+cellulose';'BSM+cellulose'};
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Experiments
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'BiomassExp'; ax = axes(fig);
semilogy(ax,exp_bsm.time,exp_bsm.biomass,'ko-', 'LineWidth',lineWidth, 'MarkerFaceColor','k'); hold(ax,'on')
semilogy(ax,exp_ynb.time,exp_ynb.biomass,'ks--', 'LineWidth',lineWidth, 'MarkerFaceColor','k'); hold(ax,'off')
xlim(ax,[0, exp_bsm.time(end)])
ylim(ax,[1, ax.YLim(2)])
ax.YTickLabel = ax.YTick;
lh = legend(medium, 'Location','Best'); lh.Box = 'Off';
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, exp_bsm.biomass_units, 'FontSize',xyLabelSize)
title(ax, microbe{1}, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\BiomassExp'],'png')

% Simulations
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'BiomassSim'; ax = axes(fig);
semilogy(ax,time,biomass{1},'k-', 'LineWidth',lineWidth, 'MarkerFaceColor','k'); hold(ax,'on')
semilogy(ax,time,biomass{2},'k--', 'LineWidth',lineWidth, 'MarkerFaceColor','k'); hold(ax,'off')
xlim(ax,[0, time(end)])
ylim(ax,[0, ax.YLim(2)])
ax.YTickLabel = ax.YTick;
lh = legend(medium, 'Location','Best'); lh.Box = 'Off';
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, microbe{1}, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\BiomassSim'],'png')






