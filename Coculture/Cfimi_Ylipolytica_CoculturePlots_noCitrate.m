
clear variables
close all
clc

%% Load Data

saveDataName = 'Cfimi_Ylipolytica_CoculturePlots_noCitrate';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% Parameters
load('parameters.mat');

% Medium
mediumName = 'YNB_noCitrate';

% C. fimi Monoculture
loadDataName = ['Cfimi_' mediumName];
mono_Cf = load([pwd '\' loadDataName '\' loadDataName '.mat']);
% Y. lipolytica Monoculture
loadDataName = ['Ylipolytica_' mediumName];
mono_Yl = load([pwd '\' loadDataName '\' loadDataName '.mat']);
% Coculture
loadDataName = ['Cfimi_Ylipolytica_' mediumName];
co_CfYl = load([pwd '\' loadDataName '\' loadDataName '.mat']);

%% Plot Variables

microbe = {'C. fimi'; 'Y. lipolytica'};
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

Cfimi_color = [84,39,143]./256;
Ylipolytica_color = [166,54,3]./256;

N = numel(co_CfYl.time);

%% Biomass Semilog Plot

n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass';
% C. fimi
ax = subplot(1,2,1); clear h
h(1) = semilogy(ax,mono_Cf.time,1E3.*mono_Cf.biomass{1},'--', 'LineWidth',lineWidth, 'Color',Cfimi_color); hold(ax, 'on');
h(2) = semilogy(ax,co_CfYl.time,1E3.*co_CfYl.biomass{1},'-', 'LineWidth',lineWidth, 'Color',Cfimi_color); hold(ax, 'off');
xlim(ax,[0, co_CfYl.time(end)])
ylim(ax,[0, 2E3])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
lh = legend(h,{'Monoculture','Coculture'}, 'Location','Best'); lh.Box = 'Off';
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, microbe{1}, 'FontSize',titleSize)
% Y. lipolytica
ax = subplot(1,2,2);
h(1) = semilogy(ax,mono_Yl.time,1E3.*mono_Yl.biomass{1},'--', 'LineWidth',lineWidth, 'Color',Ylipolytica_color); hold(ax, 'on');
h(2) = semilogy(ax,co_CfYl.time,1E3.*co_CfYl.biomass{2},'-', 'LineWidth',lineWidth, 'Color',Ylipolytica_color); hold(ax, 'off');
xlim(ax,[0, co_CfYl.time(end)])
ylim(ax,[0, 2E3])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
lh = legend(h,{'Monoculture','Coculture'}, 'Location','Best'); lh.Box = 'Off';
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, microbe{2}, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\Biomass'],'png')

%% Crossfeedograms

markerSize = 15;
arrowHead = 4;

% No Citrate
[co_CfYl.rxns,co_CfYl.idx1,co_CfYl.idx2] = intersect(co_CfYl.exchFlux_rxns{1},co_CfYl.exchFlux_rxns{2});
rm_idx = find(sum(abs(co_CfYl.exchFlux_mean{1}(:,co_CfYl.idx1)))==0 & sum(abs(co_CfYl.exchFlux_mean{2}(:,co_CfYl.idx2)))==0);
co_CfYl.rxns(rm_idx) = []; co_CfYl.idx1(rm_idx) = []; co_CfYl.idx2(rm_idx) = [];
[~,rm_idx,~] = intersect(co_CfYl.rxns,{'H+','H2O'});
co_CfYl.rxns(rm_idx) = []; co_CfYl.idx1(rm_idx) = []; co_CfYl.idx2(rm_idx) = [];
[mono_Cf.rxns,mono_Cf.idx1,~] = intersect(mono_Cf.exchFlux_rxns{1},co_CfYl.rxns);
[mono_Yl.rxns,mono_Yl.idx1,~] = intersect(mono_Yl.exchFlux_rxns{1},co_CfYl.rxns);

% Ensure That Include All Reactions of Interest
if numel(union(co_CfYl.rxns,mono_Cf.rxns)) > numel(mono_Cf.rxns)
    [rxn,idx] = setdiff(co_CfYl.rxns,mono_Cf.rxns);
    [~,~,temp] = intersect(rxn,mono_Cf.exchFlux_rxns{1});
    mono_Cf.rxns = [mono_Cf.rxns(1:(idx-1)); rxn; mono_Cf.rxns(idx:end)];
    mono_Cf.idx1 = [mono_Cf.idx1(1:(idx-1)); temp; mono_Cf.idx1(idx:end)];
end
if numel(union(co_CfYl.rxns,mono_Cf.rxns)) > numel(co_CfYl.rxns)
    [rxn,idx] = setdiff(co_CfYl.rxns,mono_Cf.rxns);
    [~,~,temp] = intersect(rxn,co_CfYl.exchFlux_rxns{1});
    co_CfYl.rxns = [co_CfYl.rxns(1:(idx-1)); rxn; co_CfYl.rxns(idx:end)];
    co_CfYl.idx1 = [co_CfYl.idx1(1:(idx-1)); temp; co_CfYl.idx1(idx:end)];
end
if numel(union(co_CfYl.rxns,mono_Yl.rxns)) > numel(mono_Yl.rxns)
    [rxn,idx] = setdiff(co_CfYl.rxns,mono_Yl.rxns);
    [~,~,temp] = intersect(rxn,mono_Yl.exchFlux_rxns{1});
    mono_Yl.rxns = [mono_Yl.rxns(1:(idx-1)); rxn; mono_Yl.rxns(idx:end)];
    mono_Yl.idx1 = [mono_Yl.idx1(1:(idx-1)); temp; mono_Yl.idx1(idx:end)];
end
if numel(union(co_CfYl.rxns,mono_Yl.rxns)) > numel(co_CfYl.rxns)
    [rxn,idx] = setdiff(co_CfYl.rxns,mono_Yl.rxns);
    [~,~,temp] = intersect(rxn,co_CfYl.exchFlux_rxns{1});
    co_CfYl.rxns = [co_CfYl.rxns(1:(idx-1)); rxn; co_CfYl.rxns(idx:end)];
    co_CfYl.idx1 = [co_CfYl.idx1(1:(idx-1)); temp; co_CfYl.idx1(idx:end)];
end

mono_Cf.fluxType = zeros(N,numel(mono_Cf.rxns));
mono_Yl.fluxType = zeros(N,numel(mono_Yl.rxns));
co_CfYl.fluxType = zeros(N,numel(co_CfYl.rxns));

% Crossfeedogram - Coculture
plotFolderName = 'Coculture';
if isdir([pwd '\' saveDataName '\' plotFolderName])==0; mkdir([pwd '\' saveDataName '\' plotFolderName]); end
n = 2; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Crossfeedogram';
for r = 1:numel(co_CfYl.rxns)
    cf_limits = 1.1.*max(max(abs(co_CfYl.exchFlux_mean{1}(:,co_CfYl.idx1(r))))).*[-1,1];
    if abs(diff(cf_limits)) <= 1E-3; cf_limits = [cf_limits(1)-1; cf_limits(2)+1]; end
    yl_limits = 1.1.*max(max(abs(co_CfYl.exchFlux_mean{2}(:,co_CfYl.idx2(r))))).*[-1,1];
    if abs(diff(yl_limits)) <= 1E-3; yl_limits = [yl_limits(1)-1; yl_limits(2)+1]; end
    ax = subplot(4,3,r);
    plot(ax,cf_limits(1):diff(cf_limits)/2:cf_limits(2),zeros(size(cf_limits(1):diff(cf_limits)/2:cf_limits(2))),'k--', 'LineWidth',0.5*lineWidth); hold(ax,'on');
    plot(ax,zeros(size(yl_limits(1):diff(yl_limits)/2:yl_limits(2))),yl_limits(1):diff(yl_limits)/2:yl_limits(2),'k--', 'LineWidth',0.5*lineWidth);
    [~,~,~,~,~,co_CfYl.fluxType(:,r)] = crossfeedogram(co_CfYl.exchFlux_mean{1}(:,co_CfYl.idx1(r)),co_CfYl.exchFlux_mean{2}(:,co_CfYl.idx2(r)),3,fig,ax,markerSize,arrowHead,0.5); hold(ax,'off');
    grid(ax,'on'); axis(ax,'square','tight'); ax.FontSize = f*axesLabelSize;
    ax.XLabel.String = []; ax.YLabel.String = [];
    title(ax,co_CfYl.rxns(r), 'FontSize',f*titleSize)
end
[~,hx] = suplabel([microbe{1} ' Exchange Flux']); set(hx, 'FontSize',xyLabelSize);
[~,hy] = suplabel([microbe{2} ' Exchange Flux'],'y'); set(hy, 'FontSize',xyLabelSize);
saveas(fig,[pwd '\' saveDataName '\' plotFolderName '\Crossfeedogram_subplot'],'png')

% Crossfeedogram - Monoculture
plotFolderName = 'Monoculture';
if isdir([pwd '\' saveDataName '\' plotFolderName])==0; mkdir([pwd '\' saveDataName '\' plotFolderName]); end
n = 3; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Crossfeedogram';
for r = 1:numel(co_CfYl.rxns)
    cf_limits = 1.1.*max(max(abs(mono_Cf.exchFlux_mean{1}(:,mono_Cf.idx1(r))))).*[-1,1];
    if abs(diff(cf_limits)) <= 1E-3; cf_limits = [cf_limits(1)-1; cf_limits(2)+1]; end
    yl_limits = 1.1.*max(max(abs(mono_Yl.exchFlux_mean{1}(:,mono_Yl.idx1(r))))).*[-1,1];
    if abs(diff(yl_limits)) <= 1E-3; yl_limits = [yl_limits(1)-1; yl_limits(2)+1]; end
    ax = subplot(4,3,r);
    plot(ax,cf_limits(1):diff(cf_limits)/2:cf_limits(2),zeros(size(cf_limits(1):diff(cf_limits)/2:cf_limits(2))),'k--', 'LineWidth',0.5*lineWidth); hold(ax,'on');
    plot(ax,zeros(size(yl_limits(1):diff(yl_limits)/2:yl_limits(2))),yl_limits(1):diff(yl_limits)/2:yl_limits(2),'k--', 'LineWidth',0.5*lineWidth);
    % C. fimi
    [~,~,l,s,a,mono_Cf.fluxType(:,r)] = crossfeedogram(mono_Cf.exchFlux_mean{1}(:,mono_Cf.idx1(r)),zeros(N,1),3,fig,ax,markerSize,arrowHead,0.5); hold(ax,'on');
    l.Color = Cfimi_color; s.MarkerFaceColor = Cfimi_color; s.Marker = 'v';
    a(isnan(a)) = []; a(a==0) = []; for ii = 1:numel(a); set(a(ii), 'Color',Cfimi_color); end
    % Y. lipolytica
    [~,~,l,s,a,mono_Yl.fluxType(:,r)] = crossfeedogram(zeros(N,1),mono_Yl.exchFlux_mean{1}(:,mono_Yl.idx1(r)),3,fig,ax,markerSize,arrowHead,0.5); hold(ax,'off');
    l.Color = Ylipolytica_color; s.MarkerFaceColor = Ylipolytica_color; s.Marker = '^';
    a(isnan(a)) = []; a(a==0) = []; for ii = 1:numel(a); set(a(ii), 'Color',Ylipolytica_color); end
    grid(ax,'on'); axis(ax,'square','tight'); ax.FontSize = f*axesLabelSize;
    ax.XLabel.String = []; ax.YLabel.String = [];
    title(ax,co_CfYl.rxns(r), 'FontSize',f*titleSize)
end
[~,hx] = suplabel([microbe{1} ' Exchange Flux']); set(hx, 'FontSize',xyLabelSize);
[~,hy] = suplabel([microbe{2} ' Exchange Flux'],'y'); set(hy, 'FontSize',xyLabelSize);
saveas(fig,[pwd '\' saveDataName '\' plotFolderName '\Crossfeedogram_subplot'],'png')

c = cbrewer('qual','Set1',numel(co_CfYl.rxns));
% Crossfeedogram - Coculture
plotFolderName = 'Coculture';
if isdir([pwd '\' saveDataName '\' plotFolderName])==0; mkdir([pwd '\' saveDataName '\' plotFolderName]); end
cf_limits = 1.1.*max(max(abs(co_CfYl.exchFlux_mean{1}(:,co_CfYl.idx1)))).*[-1,1];
if abs(diff(cf_limits)) <= 1E-3; cf_limits = [cf_limits(1)-1; cf_limits(2)+1]; end
yl_limits = 1.1.*max(max(abs(co_CfYl.exchFlux_mean{2}(:,co_CfYl.idx2)))).*[-1,1];
if abs(diff(yl_limits)) <= 1E-3; yl_limits = [yl_limits(1)-1; yl_limits(2)+1]; end
n = 4; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Crossfeedogram'; ax = axes(fig); clear h
plot(ax,cf_limits(1):diff(cf_limits)/2:cf_limits(2),zeros(size(cf_limits(1):diff(cf_limits)/2:cf_limits(2))),'k--', 'LineWidth',lineWidth); hold(ax,'on');
plot(ax,zeros(size(yl_limits(1):diff(yl_limits)/2:yl_limits(2))),yl_limits(1):diff(yl_limits)/2:yl_limits(2),'k--', 'LineWidth',lineWidth);
for r = 1:numel(co_CfYl.rxns)
    h(r) = scatter(ax,co_CfYl.exchFlux_mean{1}(1:3:end,co_CfYl.idx1(r)),co_CfYl.exchFlux_mean{2}(1:3:end,co_CfYl.idx2(r)), 'MarkerEdgeColor',c(r,:), 'MarkerFaceColor',c(r,:));
    plot(ax,co_CfYl.exchFlux_mean{1}(1:3:end,co_CfYl.idx1(r)),co_CfYl.exchFlux_mean{2}(1:3:end,co_CfYl.idx2(r)), 'LineWidth',lineWidth, 'Color',c(r,:));
end; hold(ax,'off');
grid(ax,'on'); axis(ax,'square','tight'); ax.FontSize = axesLabelSize;
lh = legend(h,co_CfYl.exchFlux_rxns{1}(co_CfYl.idx1), 'Location','Best'); lh.Box = 'Off';
xlabel([microbe{1} ' Exchange Flux'], 'FontSize',xyLabelSize);
ylabel([microbe{2} ' Exchange Flux'], 'FontSize',xyLabelSize);
saveas(fig,[pwd '\' saveDataName '\' plotFolderName '\Crossfeedogram'],'png')

% Crossfeedogram - Monoculture
plotFolderName = 'Monoculture';
if isdir([pwd '\' saveDataName '\' plotFolderName])==0; mkdir([pwd '\' saveDataName '\' plotFolderName]); end
cf_limits = 1.1.*max(max(abs(mono_Cf.exchFlux_mean{1}(:,mono_Cf.idx1)))).*[-1,1];
if abs(diff(cf_limits)) <= 1E-3; cf_limits = [cf_limits(1)-1; cf_limits(2)+1]; end
yl_limits = 1.1.*max(max(abs(mono_Yl.exchFlux_mean{1}(:,mono_Yl.idx1)))).*[-1,1];
if abs(diff(yl_limits)) <= 1E-3; yl_limits = [yl_limits(1)-1; yl_limits(2)+1]; end
n = 5; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Crossfeedogram'; ax = axes(fig); clear h
plot(ax,cf_limits(1):diff(cf_limits)/2:cf_limits(2),zeros(size(cf_limits(1):diff(cf_limits)/2:cf_limits(2))),'k--', 'LineWidth',lineWidth); hold(ax,'on');
plot(ax,zeros(size(yl_limits(1):diff(yl_limits)/2:yl_limits(2))),yl_limits(1):diff(yl_limits)/2:yl_limits(2),'k--', 'LineWidth',lineWidth);
for r = 1:numel(mono_Cf.rxns)
    h(r) = scatter(ax,mono_Cf.exchFlux_mean{1}(1:3:end,mono_Cf.idx1(r)),zeros(numel(1:3:N),1), 'MarkerEdgeColor',c(r,:), 'MarkerFaceColor',c(r,:), 'Marker','v');
    plot(ax,mono_Cf.exchFlux_mean{1}(1:3:end,mono_Cf.idx1(r)),zeros(numel(1:3:N),1), 'LineWidth',lineWidth, 'Color',c(r,:));
    g(r) = scatter(ax,zeros(numel(1:3:N),1),mono_Yl.exchFlux_mean{1}(1:3:end,mono_Yl.idx1(r)), 'MarkerEdgeColor',c(r,:), 'MarkerFaceColor',c(r,:), 'Marker','^');
    plot(ax,zeros(numel(1:3:N),1),mono_Yl.exchFlux_mean{1}(1:3:end,mono_Yl.idx1(r)), 'LineWidth',lineWidth, 'Color',c(r,:));
end; hold(ax,'off');
grid(ax,'on'); axis(ax,'square','tight'); ax.FontSize = axesLabelSize;
lh = legend(h,mono_Cf.exchFlux_rxns{1}(mono_Cf.idx1), 'Location','Best'); lh.Box = 'Off';
xlabel([microbe{1} ' Exchange Flux'], 'FontSize',xyLabelSize);
ylabel([microbe{2} ' Exchange Flux'], 'FontSize',xyLabelSize);
saveas(fig,[pwd '\' saveDataName '\' plotFolderName '\Crossfeedogram'],'png')

%% Metabolites

if isdir([pwd '\' saveDataName '\Metabolites'])==0; mkdir([pwd '\' saveDataName '\Metabolites']); end
cmap = [[233,233,233];[0,69,41];[21,122,62];[81,181,102];[8,48,107];[22,102,174];[79,156,204];[254,190,72];[208,112,12]]./256;
% Metabolite Concentration
for r = 1:numel(co_CfYl.rxns)
    [~,~,idx_co] = intersect(co_CfYl.rxns(r),co_CfYl.exchMets_names);
    [~,~,idx_cf] = intersect(co_CfYl.rxns(r),mono_Cf.exchMets_names);
    [~,~,idx_yl] = intersect(co_CfYl.rxns(r),mono_Yl.exchMets_names);
    n = 5+r;
    if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
    fig = figure(n); fig.Name = co_CfYl.rxns{r};
    ax1 = axes(fig); ax2 = axes(fig);
    linkaxes([ax1,ax2],'x'); % link x axes
    x = co_CfYl.time([1,end]);
    y = 1:6;
    C = [mono_Cf.fluxType(:,r),co_CfYl.fluxType(:,r),mono_Yl.fluxType(:,r)]';
    imagesc(ax1,x,y,C)
    colormap(ax1,cmap); caxis(ax1,[0, 8]); cbar = colorbar(ax1);
    cbar.Label.String = 'Exchange Type'; cbar.Ticks = 0.45:8/9:8;
    cbar.TickLabels = {'\rightarrow','Cf,Yl\rightarrow','Cf\rightarrow','Yl\rightarrow',...
        '\rightarrowCf,Yl','\rightarrowCf','\rightarrowYl','Yl\rightarrowCf','Cf\rightarrowYl'};
    h(1) = plot(ax2,co_CfYl.time,co_CfYl.exchMets_mM(:,idx_co),'k-', 'LineWidth',lineWidth); hold(ax2,'on');
    h(2) = plot(ax2,mono_Cf.time,mono_Cf.exchMets_mM(:,idx_cf),'--', 'LineWidth',lineWidth, 'Color',Cfimi_color);
    h(3) = plot(ax2,mono_Yl.time,mono_Yl.exchMets_mM(:,idx_yl),'--', 'LineWidth',lineWidth, 'Color',Ylipolytica_color); hold(ax2,'off');
    xlim(ax1,[0, co_CfYl.time(end)])
    ax2.Position = ax1.Position;
    grid(ax1,'on'); ax1.FontSize = axesLabelSize; ax2.FontSize = axesLabelSize;
    ax1.Color = 'none'; ax2.Color = 'none'; ax1.YTick = []; ax2.XTick = [];
    xlabel(ax1, 'Time [hours]', 'FontSize',xyLabelSize)
    ylabel(ax2, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
    title(ax1,co_CfYl.rxns(r), 'FontSize',titleSize)
    saveas(fig,[pwd '\' saveDataName '\Metabolites\' co_CfYl.rxns{r}],'png')
end


