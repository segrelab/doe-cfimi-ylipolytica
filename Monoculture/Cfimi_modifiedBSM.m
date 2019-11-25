
clear variables
close all
clc

%% Load Data

saveDataName = 'Cfimi_BSM';
if isdir([pwd '\' saveDataName])==0; mkdir([pwd '\' saveDataName]); end

% Experimental Data
exp = load([saveDataName '_exp.mat']);

% Model
load('..\Models\Cfimi_KBase_gfMPIPES_bigg_trsptrs.mat')
model{1} = Cfimi;

% Medium
load('..\Media\Media_modifiedBSM.mat')
medium = arrayfun(@(x) ['EX_' modified_BSM(x).id '_e'],1:numel(modified_BSM),'Uni',false);
media_names = arrayfun(@(x) [modified_BSM(x).id '_e'], 1:length({modified_BSM.id}), 'UniformOutput',false);
media_mM = 1E3.*[modified_BSM(:).M]; % [mM] (1 M * 1E3 mM/1 M)
% Essential Metabolites
load('..\Media\essentialMets_Cfimi.mat')
added_names = arrayfun(@(x) [Cfimi_essentialMets(x).id '_e'], 1:length({Cfimi_essentialMets.id}), 'UniformOutput',false);
[~,idx_add] = setdiff(added_names,media_names);
media_names = [media_names added_names(idx_add)];
media_mM = [media_mM, 1E3.*(5E-5.*ones(1,length(added_names(idx_add))))]; % [mM] (1 M * 1E3 mM/1 M)
% Cellulose
media_names = [media_names 'cellulose_e'];
media_mM = [media_mM, 5.6]; % [mM]

%% Set Exchange Bounds

for numModel = 1:length(model)
    % Find Exchange Reactions
    exch_rxns = findExchRxns(model{numModel});
    % Change Bounds
    model{numModel}.lb(exch_rxns) = -100;
    [~,glc_idx] = intersect(model{numModel}.rxns,'EX_glc__D_e');
    model{numModel}.lb(glc_idx) = -10;
    [~,o2_idx] = intersect(model{numModel}.rxns,'EX_o2_e');
    model{numModel}.lb(o2_idx) = -20;
end

%% Define Parameter Values

biomass_0 = 1E-5.*ones(size(model)); % gCDW
dt = 0.1; % hr
Tend = 48; % hrs
N = Tend/dt;
volume = 250E-3; % L
Km = 0.01.*ones(size(model)); % mM
Vmax = 10.*ones(size(model)); % mmol/gCDW/hr
max_biomass = 2.2; % gCDW
enzymeModel = 1;

%% Run dFBA

[time,biomass,flux,exchMets_amt,exchMets_names,feasibilityFlag] = dFBA_cellulase(model,media_names,media_mM.*volume,1,biomass_0,dt,N,volume,Km,Vmax,max_biomass,enzymeModel);
exchMets_mM = exchMets_amt./volume; % [mM] (mmol/L)
exchMets_names = strrep(strrep(erase(exchMets_names,'_e'),'__','-'),'_','-');

save([pwd '\' saveDataName '\' saveDataName '.mat'],'time','biomass','flux','exchMets_amt','exchMets_names','feasibilityFlag');

%% Plots

% Remove O2, H+, and H2O because very large values
[~,o2_idx,~] = intersect(exchMets_names,'o2');
[~,h_idx,~] = intersect(exchMets_names,'h');
[~,h2o_idx,~] = intersect(exchMets_names,'h2o');
exchMets_mM(:,[o2_idx, h_idx, h2o_idx]) = [];
exchMets_names([o2_idx, h_idx, h2o_idx]) = [];

% Enzyme
[~,cellulase_idx,~] = intersect(exchMets_names,'cellulase');
if ~isempty(cellulase_idx); cellulase_mM = exchMets_mM(:,cellulase_idx); end
exchMets_mM(:,cellulase_idx) = [];
exchMets_names(cellulase_idx) = [];

% Metabolite Indices
[~,glc_idx,~] = intersect(exchMets_names,'glc-D');
[~,cellulose_idx,~] = intersect(exchMets_names,'cellulose');
[~,n_idx,~] = intersect(exchMets_names,'no3');
[~,p_idx,~] = intersect(exchMets_names,'pi');
[~,s_idx,~] = intersect(exchMets_names,'so4');

microbe = {'C. fimi'};
medium = 'BSM+cellulose';
xyLabelSize = 30;
axesLabelSize = 24;
titleSize = 30;
lineWidth = 3;

% Biomass Semilog Plot
n = 1;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'BiomassExp'; ax = axes(fig);
semilogy(ax,exp.time,exp.biomass,'ko-', 'LineWidth',lineWidth, 'MarkerFaceColor','k');
xlim(ax,[0, exp.time(end)])
ylim(ax,[0, ax.YLim(2)])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, exp.biomass_units, 'FontSize',xyLabelSize)
title(ax, medium, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\BiomassExp'],'png')

f = 0.6;
% Biomass, Glucose, Cellulose, Nitrate, Phosphate, Sulfate, and Enzyme
n = 2;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass_Metabolites_Enzyme';
% biomass
ax = subplot(3,1,1);
semilogy(ax,time,1E3.*biomass{1},'k-', 'LineWidth',lineWidth);
xlim(ax,[0, time(end)])
ylim(ax,[0, 2.5])
ax.XTickLabel = []; ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ylabel(ax, {'Total'; 'Biomass'; '[mg CDW]'}, 'FontSize',f*xyLabelSize)
title(ax, medium, 'FontSize',titleSize)
% cellulose, glucose, n, p, s
ax = subplot(3,1,2);
plot(ax,time,exchMets_mM(:,cellulose_idx),'k-','LineWidth',lineWidth); hold(ax,'on')
plot(ax,time,exchMets_mM(:,glc_idx),'--','LineWidth',lineWidth, 'Color',0.33.*[1,1,1]);
plot(ax,time,exchMets_mM(:,n_idx),'b:','LineWidth',lineWidth);
plot(ax,time,exchMets_mM(:,p_idx),'g:','LineWidth',lineWidth);
plot(ax,time,exchMets_mM(:,s_idx),'r:','LineWidth',lineWidth); hold(ax,'off')
xlim(ax,[0, time(end)])
ax.XTickLabel = []; 
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
ylabel(ax, {'Metabolite'; 'Concentration'; '[mM]'}, 'FontSize',f*xyLabelSize)
lh = legend({'Cellulose','Glucose','Nitrate','Phosphate','Sulfate'}, 'Location','Best'); lh.Box = 'Off';
% enzyme
ax = subplot(3,1,3);
plot(ax,time,cellulase_mM,'k-','LineWidth',lineWidth)
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = f*axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',f*xyLabelSize)
ylabel(ax, {'Cellulase'; 'Concentration'; '[mM]'}, 'FontSize',f*xyLabelSize)
saveas(fig,[pwd '\' saveDataName '\Biomass_Metabolites_Enzyme'],'png')

% Biomass Semilog Plot
n = 3;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Biomass'; ax = axes(fig);
semilogy(ax,time,1E3.*biomass{1},'k-', 'LineWidth',lineWidth);
xlim(ax,[0, time(end)])
ylim(ax,[0, 2.5])
ax.YTickLabel = ax.YTick;
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Total Biomass [mg CDW]', 'FontSize',xyLabelSize)
title(ax, medium, 'FontSize',titleSize)
saveas(fig,[pwd '\' saveDataName '\Biomass'],'png')

% Metabolites
[~,col] = find(abs(diff(exchMets_mM)) > 1E-9); col = unique(col);
c = cbrewer('qual','Set1',numel(col));
n = 4;
if ishandle(n); clf(findobj('Type','Figure', 'Number',n)); end
fig = figure(n); fig.Name = 'Metabolites'; ax = axes(fig);
set(ax, 'ColorOrder',c, 'NextPlot','replacechildren');
h = plot(ax,time,exchMets_mM(:,col),'-', 'LineWidth',lineWidth);
xlim(ax,[0, time(end)])
grid(ax,'on'); ax.FontSize = axesLabelSize;
xlabel(ax, 'Time [hours]', 'FontSize',xyLabelSize)
ylabel(ax, 'Metabolite Concentration [mM]', 'FontSize',xyLabelSize)
title(ax, medium, 'FontSize',titleSize)
lh = legend(h,exchMets_names(col), 'Location','Best'); lh.Box = 'Off';
saveas(fig,[pwd '\' saveDataName '\Metabolites'],'png')





