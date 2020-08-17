%% Fig 2 shows factors that affect expansion of seizure territory
%% Simulate it
Exp2;
v = R.Var.V - R.Var.phi;
rate = O.param.f(v);
rate = rate/O.param.f_max(1);
%% Panel A
fig = figure('Units', 'inches', 'Position', [1 1 8 6]);
ax = axes('NextPlot','add','Position',[0.05 0.55 0.85 0.35],'YDir','normal','Box','off','Color','none','clipping','off');
imagesc(ax(1),[0 85],[0 1],rate(:,1:85000));
ylabel(ax(1),'Space');
ax.XLim = [0 85];
ax.YLim = [0 1];
ax.XTick = 25:25:100;
ax.YTick = [0 1];
ax.XAxis.Visible = 'off';
% Set color
H = hot(360);
colormap(ax(1),H);
cbar = colorbar(ax(1),'eastoutside');
ylabel(cbar,'Normalized firing rate');
cbar.Position = ax.Position;cbar.Position(1) = 0.91;cbar.Position(3) = 0.02;
caxis(ax,[0 0.6]);
cbar.YTick = [0 0.6];
% Marker: stimulation
Ls = plot(ax,[2 5 5 2 2],[0.1 0.1 0.15 0.15 0.1],'Color',[0 1 0]); % External stimulation
% Marker: time scale 
L_scale = plot(ax,ax.XLim(1)+[0 10],[-0.05 -0.05],'LineWidth',2,'Color',[0 0 0]);
% Marker: zoom-in
L_zoom = plot(ax,[60 60.5],[-0.05 -0.05],'Color',[0 0 0],'LineWidth',2);
% Marker: stages
Lp = plot(ax(1),ax(1).XLim,[1 1],'Color',[0    0.4471    0.7412],'LineWidth',2);
Lp(2) = scatter(ax(1),[5+5 18 67 76],[1 1 1 1],48,[0    0.4471    0.7412],'filled','Marker','diamond');
Lp(2).Marker='diamond';
TxT = text(0,1.12,{'Pre','ictal'});
TxT(2) = text(14,1.12,{'Ictal','tonic'});
TxT(3) = text(42,1.12,{'Ictal','clonic'});
TxT(4) = text(71.5,1.12,{'Pre','termination'});
TxT(5) = text(80,1.12,{'Post','ictal'});
TxT(6) = text(88,1.12,{'<stages>'});
set(TxT,'Color',[0    0.4471    0.7412],'HorizontalAlignment','center');
%% Panel B
ax(2) = axes(fig);
ax(2).Position =  [0.05 0.1 0.1 0.35];
C(2) = imagesc(ax(2),[60 60.5],[0 1],rate(:,60000:60500));
ax(2).YDir = 'normal';
ax(2).YTick = [0 1];
ax(2).XTick = [];
colormap(ax(2),H);
caxis(ax(2),[0 0.6]);
xlabel(ax(2),'0.5 sec','Color',[0.15 0.15 0.15]);
ylabel(ax(2),'Space');
set(ax,'NextPlot','add');
%% Panel C
ax(3) = axes(fig);
ax(3).Position =  [0.22 0.1 0.1 0.35];
L_Cl = plot(ax(3),R.Var.Cl_in(:,60000),linspace(0,1,max(O.n)),'LineWidth',1,'Color',[0.8500    0.3300    0.1000]);
ax(3).Color = 'none';
ax(3).XAxis.Color = [0.8500    0.3300    0.1000];
ax(3).YTick = [];
axis(ax(3),'tight');
xlabel(ax(3),'Cl_i_n (mM)');
ax(3).XTick = [8 12 16];
ax(4) = axes(fig);
ax(4).Position =  ax(3).Position;
L_K = plot(ax(4),R.Var.g_K(:,60000)/O.param.f_max(1),linspace(0,1,max(O.n)),'LineWidth',1,'Color',[0    0.4500    0.7400]);
ax(4).Color = 'none';
ax(4).Box = 'off';
ax(4).XAxis.Color = [0    0.4500    0.7400];
ax(4).XAxisLocation = 'top';
ax(4).YTick = [];
axis(ax(4),'tight');
ax(4).XLim(1) = 0;
ax(4).XTick = [0 5];
xlabel(ax(4),'g_K (nS)');
set(ax(3:4),'Box','off');
%%
savefig(fig,'D:\Paper Drafts\Model paper\Figures\2ABC.fig');
%% Panel D, robustness of inhibition determines propagation success & speed
clearvars -except ax fig
close all;
PanelFlag = 'D';
Condition = 3000:1000:8000; % tau_Cl
for n_trial = 1:numel(Condition)
    Exp2;
    % find seizure territory by convex hull
    v = R.Var.V - R.Var.phi;
    v = v(:,5001:round(O.t));    % Only count after external input has been shut off
    rate = O.param.f(v)/mean(O.param.f_max(:));
    idx = rate > 0.1;
    idx = diff([idx(1,:);idx]);
    [Y,X] = find(idx~=0);
    K = convhull(Y,X);
    Y = Y(K);X=X(K);
    sel = diff([0;Y])>0 & diff([0;X])>0;
    Y = Y(sel); % Only reach the expansion at the upper border
    X = X(sel);
    Rec_Y{n_trial} = Y/max(O.n); % Normalized space
    Rec_X{n_trial} = X/1000; % time unit: second
    close all;    
    clearvars -except n_trial Condition PanelFlag Rec_X Rec_Y
end
fig=figure;
ax = axes;ax.NextPlot = 'add';ax.Color='none';
for n_trial = 1:numel(Condition)
    l(n_trial) = plot(Rec_X{n_trial}/1000,Rec_Y{n_trial}/500,'LineWidth',1);
    l(n_trial).Color = n_trial/10 + [0.2 0.2 0.2];    
end
P = plot(ax,[0.5 0.5], [0.1,0.15],'Color',[1 0 0],'LineWidth',2);
L = legend(l,cellfun(@num2str,(num2cell(Condition/1000))'));    
L.String{1} = '3 (failed)';
L.Title.String = '\tau_C_l (second)'; 
tx = text(2,0.125,'Ext. stim.');tx.Color = [1 0 0];
L.EdgeColor = 'none';
L.Color = 'none';
ax.XLim = [0.5,100];
ax.YLim = [0 1];
ax.YTick = [0 1];
ax.XTick = 0:25:100;
xlabel(ax,'Second');
title(ax,'Expansion of seizure territory');
savefig(fig,'D:\Paper Drafts\Model paper\Figures\2D.fig');

%% Panel E, Potassium conductance controls tonic-to-clonic transition
clear;close all;
PanelFlag = 'E';
Condition = 30:2.5:50; % g_K_max
for n_trial = numel(Condition)
    Exp2;
    center_cell = 62; 
    v = R.Var.V(center_cell,:) - R.Var.phi(center_cell,:);
    v = v(:,5001:round(O.t));    % Only count after external input has been shut off
    rate = 1./(1+exp(-v/O.param.beta(1))); % Normalized firing rate
    [Pk,Idx] = findpeaks(diff([0,rate]),'MinPeakHeight',0.01);
    if isempty(Idx) % If the system stays forever in tonic stage, then just give inf
        Rec_TC_Transition(n_trial) = inf;
    else
        Idx(Idx < 1000) = []; % Remove transient
        Rec_TC_Transition(n_trial) = Idx(1);
    end
    close all;    
    clearvars -except n_trial Condition PanelFlag Rec_TC_Transition;
end
fig=figure;
ax = axes;ax.NextPlot = 'add';ax.Color='none';
Condition = Condition / 0.2 / 1000; % Change to the unit of Delta_K 
Rec_TC_Transition = Rec_TC_Transition/1000; % Change unit to second
L = plot(ax,Condition,Rec_TC_Transition,'LineWidth',1);
L.Marker = '.';
L.MarkerSize = 20;
xlabel(ax,'\Delta_K (nS/Hz)');
L(2) = plot(Condition(Rec_TC_Transition == inf),30,'.');
L(2).Color = L(1).Color;
L(2).MarkerSize = 20;
ylabel(ax,'Second');
title(ax,'Duration of tonic stage')
ax.XTick = 0.15:0.05:0.25;xlim(ax,[0.15 0.25]);
ax.YTick = 0:10:20;
Tx = text(Condition(1),30,'No transition','Color',[0.15 0.15 0.15]);
ax.YLim = [0 30];
savefig(fig,'D:\Paper Drafts\Model paper\Figures\2E.fig');

%% Figure 2 Assembly
fig = openfig('D:\Paper Drafts\Model paper\Figures\2ABC.fig');
fD = openfig('D:\Paper Drafts\Model paper\Figures\2D.fig');
fE = openfig('D:\Paper Drafts\Model paper\Figures\2E.fig');
ax_ref = fig.Children(3);
o = copyobj(fD.Children,fig);
ax = findobj(o,'Type','Axes');
ax.Position = ax_ref.Position;
ax.Position(1) = 0.38;
ax.Position(3) = 0.25;
ax.YLabel = [];
lo = findobj(o,'Type','Legend');
lo.Location = 'best';
o = copyobj(fE.Children,fig);
ax = findobj(o,'Type','Axes');
ax.Position = ax_ref.Position;
ax.Position(1) = 0.71;
ax.Position(3) = 0.25;
ax.YLabel = [];
lo = findobj(o,'Type','Legend');
lo.Location = 'best';
%%
export_fig Figure2 -png -transparent -r600 -c[0,0,0,0];