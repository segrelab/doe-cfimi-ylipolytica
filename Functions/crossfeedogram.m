function [fig_handle,axes_handle,line_handle,scatter_handle,arrow_handle,fluxType] = crossfeedogram(model1_flux,model2_flux,dt,fig_h,axes_h,markerSize,arrowHead,lineWidth)
%crossfeedogram Plots the flux of model 2 versus the flux of model 1 with
%arrows pointing in the direction of flux change from time point t to time
%point t+dt
%   [fig_handle,axes_handle,line_handle,scatter_handle,arrow_handle,fluxType] = crossfeedogram(model1_flux,model2_flux)
%   [fig_handle,axes_handle,line_handle,scatter_handle,arrow_handle,fluxType] = crossfeedogram(model1_flux,model2_flux,dt,fig_h,axes_h,markerSize,arrowHead,lineWidth)
%
% REQUIRED INPUTS
% model1_flux, model2_flux: matrix of extracellular metabolite flux over time [time x 1]
%
% OPTIONAL INPUTS
% dt: time step as the number of indices (default=1)
% fig_h: figure handle
% axes_h: axes handle
% markerSize: maximum markerSize (default=25)
% arrowHead: length and width of each arrow (default=5)
% lineWidth: width of the line (default=3)
%
% OUTPUTS
% fig_handle: figure handle
% axes_handle: axes handle
% line_handle: line handle
% scatter_handle: scatter plot handle
% arrow_handle: arrow annotations
% fluxType: Number indicates what type of flux occurs (same dimensions as
%  model1_flux and model2_flux):
%   0: not used by model 1 or model 2 (->)
%   1: produced by model 1 and model 2 (1,2->)
%   2: produced by model 1 and not used by model 2 (1->)
%   3: not used by model 1 and produced by model 2 (2->)
%   4: consumed by model 1 and model 2 (->1,2)
%   5: consumed by model 1 and not used by model 2 (->1)
%   6: not used by model 1 and consumed by model 2 (->2)
%   7: consumed by model 1 and produced by model 2 (2->1)
%   8: produced by model 1 and consumed by model 2 (1->2)
%
% Meghan Thommes 12/12/2019

%% Check Input Variables

% Make sure model1_flux and model2_flux are matrices
if ~ismatrix(model1_flux) || ~ismatrix(model2_flux)
    error('myfuns:crossfeedogram:incorrectInput', ...
        'Flux matrices must be matrices')
end
% Make sure flux matrices are the same size
if size(model1_flux) ~= size(model2_flux)
    error('myfuns:crossfeedogram:incorrectInput', ...
        'Flux matrices must be the same size')
end
% Make sure have enough inputs
if nargin < 2
    error('myfuns:crossfeedogram:NotEnoughInputs', ...
        'Not enough inputs: need model1_flux and model2_flux')
end
% Check for optional inputs
if ~exist('dt','var') || isempty(dt)
    dt = 1;
end
if ~exist('fig_h','var') || isempty(fig_h)
    fig_handle = figure;
else
    fig_handle = figure(fig_h);
end
if ~exist('axes_h','var') || isempty(axes_h)
    axes_handle = gca;
else
    axes_handle = axes_h;
end
if ~exist('markerSize','var') || isempty(markerSize)
    markerSize = 25;
end
if ~exist('arrowHead','var') || isempty(arrowHead)
    arrowHead = 5;
end
if ~exist('lineWidth','var') || isempty(lineWidth)
    lineWidth = 3;
end

%% Plot

% x-axis
x_limits = 1.1.*max(abs(model1_flux(:))).*[-1, 1];
if abs(diff(x_limits)) <= 1E-3; x_limits = [x_limits(1)-1; x_limits(2)+1]; end

% y-axis
y_limits = 1.1.*max(abs(model2_flux(:))).*[-1, 1];
if abs(diff(y_limits)) <= 1E-3; y_limits = [y_limits(1)-1; y_limits(2)+1]; end

% markers for relative concentration
c = cbrewer('seq','Reds',numel(model1_flux));
scatter_handle = scatter(axes_handle,model1_flux(1:dt:numel(model1_flux)),model2_flux(1:dt:numel(model1_flux)),markerSize,c(1:dt:numel(model1_flux),:),'filled');
hold(axes_handle,'on');

% line
line_handle = plot(axes_handle,model1_flux(1:dt:numel(model1_flux)),model2_flux(1:dt:numel(model1_flux)),'k-','LineWidth',lineWidth);

% arrows
arrowLength = 0.01;
arrow_handle = NaN(size(1:dt:numel(model1_flux)-dt));
for tt = 1:dt:numel(model1_flux)-dt % each time point
    x = model1_flux([tt, tt+dt]); % x-axis
    y = model2_flux([tt, tt+dt]); % y-axis
    
    % Arrow between time points 1 and 2
    p1 = [x(1),y(1)]; % first point
    v = [diff(x), diff(y)]; v(abs(v)<1E-6) = 0; % vector
    L = sqrt(sum(v.^2)); % vector length
    if 2*L > 0.1*max([x_limits(2),y_limits(2)])
        ah = annotation(fig_h,'arrow', 'HeadStyle','vback3', 'HeadLength',arrowHead, 'HeadWidth',arrowHead);
        ah.Parent = axes_handle; ah.LineStyle = 'none';
        ah.Position = [p1,0.5*v]; ah.Color = 'k';
        arrow_handle(tt) = ah;
    end
end
hold(axes_handle,'off');
xlabel(axes_handle,'Model 1 Exchange Flux')
ylabel(axes_handle,'Model 2 Exchange Flux')
axis(axes_handle,'tight')

%% Flux Type

tol = 1E-6;

% Determine if metabolites are produced, consumed, or not used
produced1 = find(model1_flux > tol); % produced by model 1
produced2 = find(model2_flux > tol); % produced by model 2
consumed1 = find(model1_flux < -tol); % consumed by model 1
consumed2 = find(model2_flux < -tol); % consumed by model 2
notused1 = find(model1_flux > -tol & model1_flux < tol); % not used by model 1
notused2 = find(model2_flux > -tol & model2_flux < tol); % not used by model 2

% Find the indices
idx_notused1_notused2 = intersect(notused1,notused2); % neither produced nor consumed
idx_pro1_pro2 = intersect(produced1,produced2); % produced by model 1 & model 2
idx_con1_con2 = intersect(consumed1,consumed2); % consumed by model 1 & model 2
idx_pro1_notused2 = intersect(produced1,notused2); % produced by model 1 & not used by model 2
idx_notused1_pro2 = intersect(produced2,notused1); % produced by model 2 & not used by model 1
idx_con1_notused2 = intersect(consumed1,notused2); % consumed by model 1 & not used by model 2
idx_notused1_con2 = intersect(consumed2,notused1); % consumed by model 2 & not used by model 1
idx_con1_pro2 = intersect(produced2,consumed1); % produced by model 2 & consumed by model 1
idx_pro1_con2 = intersect(produced1,consumed2); % produced by model 1 & consumed by model 2

% Classify exchange flux
fluxType = zeros(size(model1_flux));
fluxType(idx_notused1_notused2) = 0;
fluxType(idx_pro1_pro2) = 1;
fluxType(idx_pro1_notused2) = 2;
fluxType(idx_notused1_pro2) = 3;
fluxType(idx_con1_con2) = 4;
fluxType(idx_con1_notused2) = 5;
fluxType(idx_notused1_con2) = 6;
fluxType(idx_con1_pro2) = 7;
fluxType(idx_pro1_con2) = 8;

end

