
clear variables
close all
clc

%% Load Data

load('OD_data.mat')
wavelength = 300:50:1000;
N_Cfimi = sum(arrayfun(@(x) numel(x.Cfimi), OD_data)); % number of samples
N_Ylipolytica = sum(arrayfun(@(x) numel(x.Ylipolytica), OD_data)); % number of samples
N_mixture = sum(arrayfun(@(x) numel(x.mixture), OD_data)); % number of mixtures

saveDataName = 'Plots';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

%% Plot Variables

xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

Cfimi_color = [84,39,143]./256;
Ylipolytica_color = [166,54,3]./256;

%% Reformat Data

% Pre-Allocate
blank_ODavg = zeros(numel(wavelength),numel(OD_data));
Cfimi_ODavg = zeros(numel(wavelength),N_Cfimi);
Cfimi_ODstd = zeros(numel(wavelength),N_Cfimi);
Cfimi_OD600 = zeros(1,N_Cfimi);
Cfimi_name = cell(1,N_Cfimi);
Ylipolytica_ODavg = zeros(numel(wavelength),N_Ylipolytica);
Ylipolytica_ODstd = zeros(numel(wavelength),N_Ylipolytica);
Ylipolytica_OD600 = zeros(1,N_Ylipolytica);
Ylipolytica_name = cell(1,N_Ylipolytica);
mixture_ODavg = zeros(numel(wavelength),N_mixture);
mixture_ODstd = zeros(numel(wavelength),N_mixture);
mixture_OD600 = zeros(N_mixture,1);
mixture_name = cell(1,N_mixture);
mixture_alphaEst = zeros(N_mixture,1);
mixture_alpha = zeros(N_mixture,1);
mixture_ODest = zeros(numel(wavelength),N_mixture);
d = cell(1,N_mixture);
Cfimi_ODest = zeros(numel(wavelength),N_mixture);
Ylipolytica_ODest = zeros(numel(wavelength),N_mixture);
mixture_ODerror = zeros(numel(wavelength),N_mixture);

C_idx = 1; Y_idx = 1; M_idx = 1;
for n_expt = 1:numel(OD_data) % each experiment
    OD_data(n_expt).wavelength(OD_data(n_expt).wavelength == 999) = 1000;
    [~,w_idx,wAvg_idx] = intersect(OD_data(n_expt).wavelength,wavelength);
    
    blank_ODavg(wAvg_idx,n_expt) = nanmean(OD_data(n_expt).blank(w_idx,:),2);
    % C. fimi
    for n_dil = 1:numel(OD_data(n_expt).Cfimi) % each dilution
        Cfimi_ODavg(wAvg_idx,C_idx) = nanmean(OD_data(n_expt).Cfimi(n_dil).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),2);
        Cfimi_ODstd(wAvg_idx,C_idx) = nanstd(OD_data(n_expt).Cfimi(n_dil).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),0,2);
        Cfimi_OD600(C_idx) = Cfimi_ODavg(wavelength == 600,C_idx);
        Cfimi_name{C_idx} = num2str(round(100.*Cfimi_OD600(C_idx))./100);
        C_idx = C_idx +1;
    end
    % Y. lipolytica
    for n_dil = 1:numel(OD_data(n_expt).Ylipolytica) % each dilution
        Ylipolytica_ODavg(wAvg_idx,Y_idx) = nanmean(OD_data(n_expt).Ylipolytica(n_dil).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),2);
        Ylipolytica_ODstd(wAvg_idx,Y_idx) = nanstd(OD_data(n_expt).Ylipolytica(n_dil).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),0,2);
        Ylipolytica_OD600(Y_idx) = Ylipolytica_ODavg(wavelength == 600,Y_idx);
        Ylipolytica_name{Y_idx} = num2str(round(100.*Ylipolytica_OD600(Y_idx))./100);
        Y_idx = Y_idx +1;
    end
    % Mixture
    for n_mix = 1:numel(OD_data(n_expt).mixture) % each mixture
        mixture_ODavg(wAvg_idx,M_idx) = nanmean(OD_data(n_expt).mixture(n_mix).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),2);
        mixture_ODstd(wAvg_idx,M_idx) = nanstd(OD_data(n_expt).mixture(n_mix).OD(w_idx,:) - blank_ODavg(wAvg_idx,n_expt),0,2);
        mixture_OD600(M_idx) = mixture_ODavg(wavelength == 600,M_idx);
        mixture_alpha(M_idx) = OD_data(n_expt).mixture(n_mix).Cfimi;
        mixture_name{M_idx} = [num2str(100*mixture_alpha(M_idx)) '% C. fimi'];
        [mixture_alphaEst(M_idx),mixture_ODest(:,M_idx),Cfimi_ODest(:,M_idx),Ylipolytica_ODest(:,M_idx),d{M_idx}] = ...
            estimateOD(Cfimi_ODavg,Ylipolytica_ODavg,mixture_ODavg(:,M_idx));
        M_idx = M_idx +1;
    end
end
[~,Cfimi_order] = sort(Cfimi_OD600,'descend');
[~,Ylipolytica_order] = sort(Ylipolytica_OD600,'descend');
[~,mixture_order] = sort(mixture_alpha,'descend');
% C. fimi
Cfimi_OD600 = Cfimi_OD600(Cfimi_order);
Cfimi_ODavg = Cfimi_ODavg(:,Cfimi_order);
Cfimi_ODstd = Cfimi_ODstd(:,Cfimi_order);
Cfimi_name = Cfimi_name(Cfimi_order);
Cfimi_color = flip(cbrewer('seq','Purples',N_Cfimi+1)); Cfimi_color(end,:) = [];
% Y. lipolytica
Ylipolytica_OD600 = Ylipolytica_OD600(Ylipolytica_order);
Ylipolytica_ODavg = Ylipolytica_ODavg(:,Ylipolytica_order);
Ylipolytica_ODstd = Ylipolytica_ODstd(:,Ylipolytica_order);
Ylipolytica_name = Ylipolytica_name(Ylipolytica_order);
Ylipolytica_color = flip(cbrewer('seq','Oranges',N_Ylipolytica+1)); Ylipolytica_color(end,:) = [];
% Mixture
mixture_OD600 = mixture_OD600(mixture_order);
mixture_ODavg = mixture_ODavg(:,mixture_order);
mixture_ODstd = mixture_ODstd(:,mixture_order);
mixture_name = mixture_name(mixture_order);
mixture_alpha = mixture_alpha(mixture_order);
mixture_alphaEst = mixture_alphaEst(mixture_order);
mixture_ODest = mixture_ODest(:,mixture_order);
Cfimi_ODest = Cfimi_ODest(:,mixture_order);
Ylipolytica_ODest = Ylipolytica_ODest(:,mixture_order);
d = d(mixture_order);
mixture_cmap = cbrewer('div','PuOr',107); mixture_cmap(50:56,:) = [];
idx = round(mixture_alpha.*size(mixture_cmap,1))+1;
mixture_color = mixture_cmap(idx,:);

%% Plots

% Plot Limits
wavelength_lim = [290 1010];
Cfimi_lim = ceil(10*max(max(Cfimi_ODavg+Cfimi_ODstd)))/10;
Ylipolytica_lim = ceil(10*max(max(Ylipolytica_ODavg+Ylipolytica_ODstd)))/10;
mixture_lim = ceil(10*max(max(mixture_ODavg+mixture_ODstd)))/10;
OD_lim = [0, max([Cfimi_lim,Ylipolytica_lim])];

% OD: Monocultures
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Monoculture OD';
% C. fimi
ax = subplot(1,2,1); clear h
h = errorbar(ax,repmat(wavelength',1,N_Cfimi),Cfimi_ODavg,Cfimi_ODstd,'v-', 'LineWidth',lineWidth);
for n_dil = 1:numel(h)
    h(n_dil).Color = Cfimi_color(n_dil,:);
end
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlim(ax,wavelength_lim); ylim(ax,OD_lim)
xlabel(ax,'Wavelength [nm]', 'FontSize',xyLabelSize)
ylabel(ax,'OD', 'FontSize',xyLabelSize)
title(ax,'C. fimi', 'FontSize',titleSize)
% Y. lipolytica
ax = subplot(1,2,2); clear h
h = errorbar(ax,repmat(wavelength',1,N_Ylipolytica),Ylipolytica_ODavg,Ylipolytica_ODstd,'^-', 'LineWidth',lineWidth);
for n_dil = 1:numel(h)
    h(n_dil).Color = Ylipolytica_color(n_dil,:);
end
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlim(ax,wavelength_lim); ylim(ax,OD_lim)
xlabel(ax,'Wavelength [nm]', 'FontSize',xyLabelSize)
ylabel(ax,'OD', 'FontSize',xyLabelSize)
title(ax,'Y. lipolytica', 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\OD_Monoculture'],'png')

% OD: Mixture
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Mixture OD'; ax = axes(fig); clear h
h = errorbar(ax,repmat(wavelength',1,N_mixture),mixture_ODavg,mixture_ODstd,'o-', 'LineWidth',lineWidth);
for n_dil = 1:numel(h)
    h(n_dil).Color = mixture_color(n_dil,:);
end
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlim(ax,wavelength_lim); ylim(ax,OD_lim)
cbar = colorbar(ax); colormap(ax,mixture_cmap);
cbar.Ticks = 0:0.1:1; cbar.TickLabels = 0:10:100; cbar.Label.String = '% C. fimi';
xlabel(ax,'Wavelength [nm]', 'FontSize',xyLabelSize)
ylabel(ax,'OD', 'FontSize',xyLabelSize)
title(ax,'Mixture', 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\OD_Mixture'],'png')

% Alpha
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'alpha'; ax = axes(fig); clear h
% Distance
ax = subplot(1,2,1);
scatter(ax,mixture_alpha,arrayfun(@(x) nanmin(d{x}(:)),1:numel(d)),36,mixture_color,'filled', 'MarkerEdgeColor','k', 'LineWidth',1);
axis(ax,'square'); grid(ax,'on'); ax.FontSize = axesLabelSize;
xlim(ax,[0,1]);
xlabel(ax,'\alpha', 'FontSize',xyLabelSize)
ylabel(ax,'Distance', 'FontSize',xyLabelSize)
% Alpha
ax = subplot(1,2,2);
plot(ax,0:1,0:1,'k--', 'LineWidth',lineWidth); hold(ax,'on');
scatter(ax,mixture_alpha,mixture_alphaEst,36,mixture_color,'filled', 'MarkerEdgeColor','k', 'LineWidth',1); hold(ax,'off');
axis(ax,'square'); grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax,'\alpha', 'FontSize',xyLabelSize)
ylabel(ax,'\alpha_{est}', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Alpha'],'png')

% Residuals
n = 4;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Residuals'; ax = axes(fig); clear h
h = plot(ax,wavelength,mixture_ODest - mixture_ODavg,'o-', 'LineWidth',lineWidth);
for n_dil = 1:numel(h)
    h(n_dil).Color = mixture_color(n_dil,:);
    h(n_dil).MarkerFaceColor = mixture_color(n_dil,:);
    h(n_dil).MarkerEdgeColor = mixture_color(n_dil,:);
end
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlim(ax,wavelength_lim);
xlabel(ax,'Wavelength [nm]', 'FontSize',xyLabelSize)
ylabel(ax,'OD_{est} - OD', 'FontSize',xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Residuals'],'png')
%%
% Estimated OD
n = 5; f = 0.6;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Estimated OD'; 
for n_mix = 1:N_mixture
    ax = subplot(4,3,n_mix); clear h
    h(1) = errorbar(ax,wavelength',mixture_ODavg(:,n_mix),mixture_ODstd(:,n_mix),'ko-', 'MarkerFaceColor','k', 'LineWidth',lineWidth); hold(ax,'on')
    h(2) = plot(wavelength',mixture_ODest(:,n_mix),'s--', 'Color',mixture_color(n_mix,:), 'MarkerFaceColor',mixture_color(n_mix,:), 'LineWidth',lineWidth);
    h(3) = plot(wavelength',Cfimi_ODest(:,n_mix),'v-', 'Color',Cfimi_color(1,:), 'MarkerFaceColor',Cfimi_color(1,:), 'LineWidth',0.5*lineWidth);
    h(4) = plot(wavelength',Ylipolytica_ODest(:,n_mix),'^--', 'Color',Ylipolytica_color(1,:), 'MarkerFaceColor',Ylipolytica_color(1,:), 'LineWidth',0.5*lineWidth);
    hold(ax,'off'); grid(ax,'on'); ax.FontSize = f*axesLabelSize;
    xlim(ax,wavelength_lim); ylim(ax,OD_lim)
    title(['\alpha = ' num2str(100*mixture_alpha(n_mix)) '% & \alpha_{est} = ' num2str(100*mixture_alphaEst(n_mix)) '%'], 'FontSize',f*titleSize);
end
[~,hx] = suplabel('Wavelength [nm]'); hx.FontSize = xyLabelSize;
[~,hy] = suplabel('OD','y'); hy.FontSize = xyLabelSize;
saveas(fig,[pwd '\' saveDataName '\OD_Estimate'],'png')
%%
for n_mix = 1:N_mixture
    n = 5 + n_mix;
    if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
    fig = figure(n); fig.Name = ['Estimate' num2str(n_mix)];
    % Distance Heat Map
    ax = subplot(1,2,1); clear h
    a = linspace(0,100,size(d{n_mix},3));
    d_plot = reshape(d{n_mix},N_Cfimi*N_Ylipolytica,numel(a));
    idx = 1:N_Cfimi*N_Ylipolytica;
    [n_dil,i] = ind2sub(size(d_plot),find(d_plot == min(min(d_plot))));
    % heat map
    imagesc(ax,a,idx,d_plot)
    % highlight minimum distance
    hold(ax,'on')
    plot(ax,a(i) + [-1,1].*0.5.*unique(diff(a)),[idx(n_dil) idx(n_dil)]-0.5,'k-', 'LineWidth',lineWidth);
    plot(ax,a(i) + [-1,1].*0.5.*unique(diff(a)),[idx(n_dil+1) idx(n_dil+1)]-0.5,'k-', 'LineWidth',lineWidth);
    plot(ax,a(i) + [1,1].*0.5.*unique(diff(a)),[idx(n_dil:(n_dil+1))]-0.5,'k-', 'LineWidth',lineWidth);
    plot(ax,a(i-1) + [1,1].*0.5.*unique(diff(a)),[idx(n_dil:(n_dil+1))]-0.5,'k-', 'LineWidth',lineWidth);
    hold(ax,'off')
    colormap(ax,cbrewer('seq','Reds',256));
    axis(ax,'tight'); ax.YTick = [];
    grid(ax,'on'); ax.FontSize = axesLabelSize;
    xlabel(ax,'\alpha (% C. fimi)', 'FontSize',xyLabelSize)
    ylabel(ax,'Sample', 'FontSize',xyLabelSize)
    zlabel(ax,'Distance', 'FontSize',xyLabelSize)
    caxis(ax,[0 ceil(10*max(d_plot(:)))/10])
    cbar = colorbar(ax); cbar.Label.String = 'Distance';
    cbar.FontSize = axesLabelSize; cbar.Label.FontSize = xyLabelSize;
    % OD
    ax = subplot(1,2,2); clear h
    h(1) = errorbar(ax,wavelength',mixture_ODavg(:,n_mix),mixture_ODstd(:,n_mix),'ko-', 'MarkerFaceColor','k', 'LineWidth',lineWidth); hold(ax,'on')
    h(2) = plot(wavelength',mixture_ODest(:,n_mix),'s--', 'Color',mixture_color(n_mix,:), 'MarkerFaceColor',mixture_color(n_mix,:), 'LineWidth',lineWidth);
    h(3) = plot(wavelength',Cfimi_ODest(:,n_mix),'v-', 'Color',Cfimi_color(1,:), 'MarkerFaceColor',Cfimi_color(1,:), 'LineWidth',0.5*lineWidth);
    h(4) = plot(wavelength',Ylipolytica_ODest(:,n_mix),'^--', 'Color',Ylipolytica_color(1,:), 'MarkerFaceColor',Ylipolytica_color(1,:), 'LineWidth',0.5*lineWidth);
    hold(ax,'off');
    lh = legend(h,{'Measured Mixture','Predicted Mixture','Predicted C. fimi','Predicted Y. lipolytica'}, 'Location','NorthEast'); lh.Box = 'Off';
    xlim(ax,wavelength_lim); ylim(ax,OD_lim)
    xlabel(ax,'Wavelength [nm]', 'FontSize',xyLabelSize)
    ylabel(ax,'OD', 'FontSize',xyLabelSize)
    [~,h] = suplabel(['\alpha = ' num2str(100*mixture_alpha(n_mix)) '% & \alpha_{est} = ' num2str(100*mixture_alphaEst(n_mix)) '%'],'t');
    h.FontSize = titleSize;
    saveas(fig,[pwd '\' saveDataName '\Estimate_' num2str(n_mix)],'png')
end




