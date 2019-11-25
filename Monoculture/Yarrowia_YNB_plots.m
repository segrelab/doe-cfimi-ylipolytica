
clear variables
close all
clc

%% Load Data

saveDataName = 'Yarrowia_YNB';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% Simulations
sugar_names = {'glc__D';'arab';'man';'gal';'xyl__D'}; % glucose, arabinose, mannose, galactose, xylose
for s = 1:numel(sugar_names)
    data = load([saveDataName '\' saveDataName '_' sugar_names{s} '.mat'],'time','biomass');
    time = data.time;
    biomass(:,s) = data.biomass{1};
end
data = load([saveDataName '\' saveDataName '_all.mat'],'time','biomass');
biomass_all = data.biomass{1};

%% Plots

microbe = {'Y. lipolytica'};
medium = {'YNB+glucose';'YNB+arabinose';'YNB+mannose';'YNB+galactose';'YNB+xylose'};
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Biomass Semilog Plot - Individual Sugars
c = cbrewer('qual','Set1',numel(sugar_names));
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass_Individual'; ax = axes(fig);
set(ax, 'ColorOrder',c, 'NextPlot','replacechildren');
semilogy(ax,time,biomass,'-', 'LineWidth',lineWidth);
xlim(ax,[0, time(end)])
ylim(ax,[0, ax.YLim(2)])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, microbe{1}, 'FontSize',titleSize)
lh = legend(medium, 'Location','Best'); lh.Box = 'Off';
fig.Renderer = 'Painters';
saveas(fig,[pwd '\' saveDataName '\Biomass_Individual'],'png')

% Biomass Semilog Plot - Individual Sugars
c = cbrewer('qual','Set1',numel(sugar_names));
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass_Individual'; ax = axes(fig);
set(ax, 'ColorOrder',c, 'NextPlot','replacechildren');
semilogy(ax,time,biomass_all,'k-', 'LineWidth',lineWidth);
xlim(ax,[0, time(end)])
ylim(ax,[0, ax.YLim(2)])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, microbe{1}, 'FontSize',titleSize)
fig.Renderer = 'Painters';
saveas(fig,[pwd '\' saveDataName '\Biomass_All'],'png')






