
clear variables
close all
clc

%% Load Data - EC 3.2.1.4 cellulase

load('cellulase_values.mat');

saveDataName = 'Kinetic_Parameters';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

%% Plots

Cfimi_color = [84,39,143]./256;
Ylipolytica_color = [166,54,3]./256;
enzyme_color = [0,109,44]./256;

xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Km
Km_idx = find(arrayfun(@(x) ~isempty(intersect(Km{x,3},'Cellulomonas fimi','rows')),1:size(Km,1)));
n = 1; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Km';
% All
ax = subplot(3,1,1);
h = histogram([Km{:,1}], 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([Km{Km_idx,1}],h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = h.BinEdges;
% [0,20]
ax = subplot(3,1,2);
h = histogram([Km{:,1}], 0:20, 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([Km{Km_idx,1}], h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = 0:20;
ylabel(ax, 'Number of Entries', 'FontSize',xyLabelSize)
% [0,2]
ax = subplot(3,1,3);
h = histogram([Km{:,1}], 0:0.1:2, 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([Km{Km_idx,1}], h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = h.BinEdges;
xlabel(ax, 'K_M [mM]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Km'],'png')

% kcat
kcat_idx = find(arrayfun(@(x) ~isempty(intersect(kcat{x,3},'Cellulomonas fimi','rows')),1:size(kcat,1)));
n = 2; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'kcat';
% All
ax = subplot(3,1,1);
h = histogram([kcat{:,1}],0:100:3300, 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([kcat{kcat_idx,1}],h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = h.BinEdges;
% [0,250]
ax = subplot(3,1,2);
h = histogram([kcat{:,1}], 0:10:250, 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([kcat{kcat_idx,1}], h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = h.BinEdges;
ylabel(ax, 'Number of Entries', 'FontSize',xyLabelSize)
% [0,10]
ax = subplot(3,1,3);
h = histogram([kcat{:,1}], 0:0.25:10, 'EdgeColor','k', 'FaceColor','k', 'FaceAlpha',1); hold(ax,'on');
histogram([kcat{kcat_idx,1}], h.BinEdges, 'EdgeColor','k', 'FaceColor',Cfimi_color, 'FaceAlpha',1); hold(ax,'off');
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ax.XMinorGrid = 'on'; ax.XAxis.MinorTickValues = h.BinEdges;
xlabel(ax, 'k_{cat} [1/s]', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\kcat'],'png')

fprintf('Km average = %f mM\n', nanmean([Km{:,1}]))
fprintf('Km C. fimi average = %f mM\n', nanmean([Km{Km_idx,1}]))
fprintf('kcat average = %f 1/s or %f 1/hr\n', nanmean([kcat{:,1}]),60*nanmean([kcat{:,1}]))
fprintf('kcat C. fimi average = %f 1/s or %f 1/hr\n', nanmean([kcat{kcat_idx,1}]),60*nanmean([kcat{kcat_idx,1}]))


