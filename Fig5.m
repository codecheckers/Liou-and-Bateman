%% Prepare figure
fig = figure('Units','inches','Position',[1 1 8 5.5]);
ax = axes('Position',[0.05 0.73 0.25 0.2],'Color','none','Box','off','NextPlot','add'); % Original raster
ax.XAxis.Visible = 'off';
ax.YLim = [0.25 0.75];
ax.YTick = [0.25 0.75];

ax_lfp = copyobj(ax,fig);
ax_lfp.Position([2,4])=[0.6 0.11];

ax_W = copyobj(ax(1),fig);
ax_W.Position = ax(1).Position;
ax_W.Position(1) = 0.34;
ax_W.Position(3) = ax_W.Position(3) * fig.Position(4) / fig.Position(3);

ax_schematic = copyobj(ax(1),fig);
ax_schematic.Position(1) = ax_W.Position(1) + ax_W.Position(3) + 0.02;
ax_schematic.Position(3) = 0.08;
ax_schematic.YTick = [];

ax(2) = copyobj(ax(1),fig); % Learned raster
ax(2).Position(1) = 0.62;
ax(2).Position(3) = 0.33;
ax(2).YTick = [];

ax_lfp(2) = copyobj(ax_lfp,fig);
ax_lfp(2).Position([1 3])=ax(2).Position([1 3]);
set(ax_lfp,'YTick',[]);

% Patient LFP
ax(3) = copyobj(ax(1),fig);
ax(3).Position(2) = 0.33;
ax(3).YTick = [];

% For image frames of backward traveling waves
ax(4) = copyobj(ax(3),fig);
ax(4).Position(1) = ax_W.Position(1);
ax(4).Position(3) = 0.6; 

% For image frames of slow propagation
ax(5) = copyobj(ax(3),fig);
ax(5).Position(4) = ax(5).Position(4)/3.3;
ax(5).Position(3) = ax(5).Position(3)/2.2;
ax(5).Position(1) = ax(3).Position(1) + ax(3).Position(3) - ax(5).Position(3);

% For repetitive c5 seizures
ax(6) = copyobj(ax(1),fig);
ax(6).Position([2,3,4]) = [0.08 0.6 0.23];
ax(6).XAxis.Visible = 'on';

% For c5 invasion sequence
ax(7) = copyobj(ax(6),fig);
ax(7).Position([1,3]) = [0.72 0.23];
ax(7).YTick = [];
ax(7).XTick = [];

ylabel(ax(1),'Space');
ylabel(ax_lfp(1),'LFP');

ylabel(ax(3),'Patient LFP');
title(ax(1),'Before learning');
title(ax(2),'After learning');

title(ax(7),'Sequence comparision');
xlabel(ax(7),'Sequence');
ylabel(ax(7),'Sequence');

%% Panel A - before learning
n_trial = 1;
Exp5;
TransferSpikeRecord(R);
%% Plot Panel A
plot(R);
l_raster = R.Graph.Raster.Children.Children;
l_raster.YData = l_raster.YData/max(O.n);
l_raster.XData = l_raster.XData/1000;
l_raster.Parent = ax(1);
pa_t = repmat(stim_t/1000,[2 1]);
pa_x = [stim_x,fliplr(stim_x)];
pa = patch(ax(1),pa_t(:),pa_x,[1 0 0],'FaceAlpha',0.7,'EdgeColor','none');
t_clip = [2-eps 12];
index_clip = round(t_clip*1000);
AxesClip(ax(1),t_clip,0);
ax(1).Children = ax(1).Children(end:-1:1);
ax(1).XLim = [0 diff(t_clip)+eps];
lfp = LFP(R);
elec = [0.6];
ax_lfp(1).ColorOrderIndex = 4;
for i = 1:numel(elec)
    l_lfp(i) = plot(ax_lfp(1),(index_clip(1):index_clip(2))/1000, ...
                              -lfp(round(elec(i)*max(O.n)),index_clip(1):index_clip(2)));
    l_lfp_mark(i) = plot(ax(1),ax(1).XLim, ...
                               [elec(i),elec(i)], ...
                               '--','Color',l_lfp(i).Color,'LineWidth',1);
end
if numel(l_lfp)>1
    l_lfp(2).YData = l_lfp(2).YData + max(l_lfp(1).YData) - min(l_lfp(2).YData);
end
axis(ax_lfp,'tight');
xlim(ax_W,[0.25 0.75]);
ylim(ax_W,[0.25 0.75]);

%% Learn
KernelToMultiplication(P_E);
BatchSTDPLearning(P_E);

%% Panel B, the learned matrix
WMatrix{1} = P_E.STDP.W;
dW = P_E.STDP.W - 1;
%%
dWmax = max(abs(dW(:)));
eta_adjustment_factor = 0.2 / dWmax; % Let's force it approximately 20% by the first learning
eta_adjustment_factor = round(eta_adjustment_factor*100)/100;
dW = dW * eta_adjustment_factor;
Im_W = imagesc(ax_W,[0 1],[0 1],dW*eta_adjustment_factor);
ax_W.YDir = 'normal';
ax_W.XDir = 'reverse';
ax_W.Box = 'off';
ax_W.XTick = [];
ax_W.YTick = [];
ax_W.XAxis.Visible = 'on';
ylabel(ax_W,'Post-synaptic');
xlabel(ax_W,'Pre-synaptic');
ax_W.XAxisLocation = 'top';
l_diag = plot(ax_W,[0 1],[0 1],'Color',[0.15 0.15 0.15]);
%%
colormap(ax_W,jet(360));
cbar_W = colorbar('peer',ax_W,'South');
cbar_W.Position = ax_W.Position;
cbar_W.Position(2) = cbar_W.Position(2) - 0.05;
cbar_W.Position(4) = 0.008;
xlabel(cbar_W,'\DeltaW');


%% Panel C Re-set the field back to O.t = 1 & simulate it again 

O.t = 1;
O.V = R.Var.V(:,1);
O.phi = R.Var.phi(:,1);
O.Cl_in = R.Var.Cl_in(:,1);
O.g_K = R.Var.g_K(:,1);
O.Input.E = 0;
O.Input.I = 0;
Refresh(R);
close(R.Graph.Raster);
close(R.Graph.Figure);
% dW = P_E.STDP.W - 1;
%%
n_trial = 2;
P_E.W = P_E.W .* (1 + dW);
Exp5;
%% learning rate
TransferSpikeRecord(R);
plot(R);
l_raster(2) = R.Graph.Raster.Children.Children;
l_raster(2).YData = l_raster(2).YData/max(O.n);
l_raster(2).XData = l_raster(2).XData/1000;
l_raster(2).Parent = ax(2);
t_clip = [5,18];
AxesClip(ax(2),t_clip,0);
ax(2).XLim = [0 diff(t_clip)];
ax(2).Position(3) = ax(1).Position(3) * diff(ax(2).XLim) / diff(ax(1).XLim);
ax_lfp(2).Position(3) = ax(2).Position(3);
lfp = LFP(R);
elec = 0.6;
l_lfp(2) = plot(ax_lfp(2),((t_clip(1)*1000):(t_clip(2)*1000))/1000, ...
                          -lfp(round(elec*max(O.n)),(t_clip(1)*1000):(t_clip(2)*1000)),'Color',[0.47 0.67 0.19]);
l_lfp_mark(2) = plot(ax(2),ax(2).XLim,[elec,elec],'--','Color',l_lfp(2).Color,'LineWidth',1);
axis(ax_lfp,'tight');
xlim(ax_W,[0.25 0.75]);
ylim(ax_W,[0.25 0.75]);
%% Put in schematics
Im_schematic = imread('D:\Paper Drafts\Model paper\Figures\4Schematic.png');
image(ax_schematic,Im_schematic);
axis(ax_schematic,'image');
ax_schematic.YAxis.Color = [1 1 1];
%% Put patient data - c7 - herald spikes (this requires you to connect with the patient data mobile disk)
map = load('PatientElectrodeMap'); % The standard Utah array map for patient c1-c7.
map = map.map; 
map = map(4:end,1:9); % Tailor it to produce patient c7's map (remove damaged electrodes) 
map(1,5:end) = -1;
map(2,end) = -1; 
ch = map(map > 0);
m = matfile('E:\data\Patient Data\C7\SeizureDownSampled30Kto2KC7.mat');
LFP_c7 = m.seizure_downsampled;
%%
LFP_c7_avg = mean(LFP_c7(ch,:));
[b,a] = butter(4,50/1000); % Prevent aliasing while plotting & as defined as LFP
LFP_c7_avg = filtfilt(b,a,LFP_c7_avg);
%
T_sel = [159,168];
LFP_c7_avg = LFP_c7_avg((T_sel(1)*2000):(T_sel(2)*2000));
Tvec = T_sel(1) + (1:numel(LFP_c7_avg))/2000;
l_lfp_pt = plot(ax(3),Tvec,LFP_c7_avg,'Color',[0.47 0.67 0.19]);
axis(ax(3),'tight');
l_scalebar = plot(ax(3),T_sel(1)+[4 4 5],(ax(3).YLim(2)-eps)*[1 1 1] - [500 0 0],'LineWidth',2,'Color',[0.15 0.15 0.15]);
ax(3).ColorOrderIndex = 1;
%%
Ev = EventData.load('E:\data\Patient Data\C7\SeizureMUADataC7.mat');
%%
T_burst = [162.89, 160.9 160.608 160.3]; % This is for all herald spikes
burst_sel = [1 1 0 0]; % In the main figure, we only show 2 herald spikes
T_burst = T_burst(burst_sel); % So select these events

dT_burst = -0.025:0.005:0.035; % Reduce temporal resolution
for i = 1:numel(T_burst)
    TS_fast = Ev2TS(Ev,T_burst(i) + dT_burst,'sigma',0.01);
    ImData = InterchangeMatrixArray(TS_fast.Data','map',map);
    for j = 1:numel(dT_burst)
        ImHS(i,j) = imagesc(ax(4),[0 0.9]*0.9 +(j-1),[0 0.7]*0.9 +(i-1)*0.9,ImData(:,:,j));
    end
    [~,idx] = min(abs((l_lfp_pt.XData - T_burst(i))));
    [~,idx_sub] = min(l_lfp_pt.YData((idx-500):(idx+500)));
    idx = idx - 501 + idx_sub;
    sc_t_mark(i) = scatter(ax(3),l_lfp_pt.XData(idx), l_lfp_pt.YData(idx), 36, '^', 'filled');
    sc_Im_mark(i) = scatter(ax(4),-0.2, mean(ImHS(i,1).YData) , 36 ,sc_t_mark(i).CData, '>', 'filled');
end

%
T_grid = [164.5 167.5];
TS_slow = Ev2TS(Ev,T_grid,'sigma',1);
ImData = InterchangeMatrixArray(TS_slow.Data','map',map);
clear ImSlow
for j = 1:numel(T_grid)
    ImSlow(i,j) = imagesc(ax(5),[0 0.9]*0.9 +(j-1),[0 0.7]*0.9,ImData(:,:,j));
end
axis(ax(4:5),'tight');
l_tbar = plot(ax(3),repmat(T_grid,[2 1]),repmat(ax(3).YLim,[2 1])','LineWidth',1,'LineStyle','--','Color',[0.5 0.5 0.5]);


%%
cbar = colorbar('peer',ax(4),'West');
cbar.Ticks = [0 600];
cbar.Position = ax(4).Position;
cbar.Position([1 3]) = [0.955 0.008];
cbar.TickLabels{2} = '0.6';
ylabel(cbar(1),'spks/ms','Position',[0.92 317.8333 0]);



%% Put patient data - c5 - repeated seizures (this requires you to connect with the patient data mobile disk)
map = load('PatientElectrodeMap'); % The standard Utah array map for patient c1-c7.
map = map.map; 
ch = [     1     2     3     4     6     7     8     9    10    11    13    16    17    20    22    23    24 ...
    25    27    28    29    30    31    32    33    35    36    37    39    41    45    49    53    54 ...
    55    57    59    60    61    62    63    64    65    66    67    70    71    73    75    77    79 ...
    80    83    84    86    88    89    90    91    93    94    96]; % This is based on Cathy & Garret's UtahKey file
map(~ismember(map,ch)) = -1;
%%
for i = 1:3
    Ev = EventData.load(['E:\data\Patient Data\C5\SeizureMUADataC50' num2str(i) '.mat']);
    TS = Ev2TS(Ev,120:0.1:160,'sigma',1);
    TS.Data = TS.Data(:,ch);
    % Sort according to invasion sequence
    [~, Idx_max] = max(TS.Data);
    [~, Idx_seq{i}] = sort(Idx_max);
    % While plotting, always plot as the sequence of the first seizure
    Im(i) = imagesc(ax(6),[0 40] + 45*(i-1),[0 1],TS.Data(:,Idx_seq{1})');
end
axis(ax(6),'tight');
colormap(ax(6),hot(360));
caxis(ax(6),[0 80]);
ax(6).XTick = [];
ax(6).XTick = [];
ylabel(ax(6),'electrode');
xlabel(ax(6), 'Time (40 sec)');
ax(6).XAxis.Color = [1 1 1];
ax(6).XLabel.Color = [0.15 0.15 0.15];
ax(6).XLabel.Position = [20   -0.055   -1.0000];
for i = 1:3
    Tx(i) = text(ax(6),20 + (i-1)*45,1.1,['sz. #' num2str(i)]);
    Tx(i).HorizontalAlignment = 'center';
end
%%
cbar(2) = colorbar('peer',ax(6),'South');
cbar(2).Ticks = [0 80];
cbar(2).TickLabels{2} = '          80 spikes/sec';
cbar(2).AxisLocation = 'out';
cbar(2).Position = [0.265    0.05    0.17    0.012];
%% Invasion sequence scatter plot
ax(7).ColorOrderIndex = 2;
pair_idx = [1 2;
            1 3;
            2 3];

for i = 1:size(pair_idx,1)
    sz_X = pair_idx(i,1);    
    sz_Y = pair_idx(i,2);
    [~,I0] = ismember(ch(Idx_seq{sz_X}),ch(Idx_seq{sz_Y}));    
    sc_array(i) = scatter(ax(7),1:numel(ch),I0,6,'filled');hold on;
    % Statistics
    [ax(7).UserData.rho(i),  ax(7).UserData.pValue(i)] = corr((1:numel(ch))',I0(:),'type','spearman');
end
l_diag = plot(ax(7),[0 numel(ch)],[0 numel(ch)],'Color',[0.5 0.5 0.5]);
axis(ax(7),'tight');
leg = legend(sc_array,{'sz #2 vs #1', 'sz #3 vs #1', 'sz #3 vs #2'},'Color','none','EdgeColor','none');
leg.Title.String = ['n = ' num2str(numel(ch))];




%%

%% Add panel labels
Char = 'ABCDEFG';
for i = 1:numel(Char)
    annotation('TextBox','String',Char(i),'FontSize',14,'FontWeight','bold','EdgeColor','none')
end

%% Re-catch the figure handle
fig = openfig('D:\Paper Drafts\Model paper\Figures\5.fig');
ax = fig.Children([end,end-2,3,2]);set(ax,'NextPlot','add');
ax_lfp = fig.Children([end-1,end-3]);set(ax_lfp,'NextPlot','add');
cbar = fig.Children(1);
ax_W = fig.Children(5);set(ax_W,'NextPlot','add');
cbar_W = fig.Children(4);

%% Continuous epileptogenesis
O = SpikingModel('Exp5Template');
[ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O );
KernelToMultiplication(P_E);
P_E.STDP.Wmax = 1.5;
W0 = P_E.W;
load('dW.mat');
dWmax = 0.1; % How much you learn from the seizure exposure
eta = dWmax / max(abs(dW(:)));
P_E.STDP.W = P_E.STDP.W + dW*eta;
P_E.W = W0 .* P_E.STDP.W;
%% Every 1 minute you learn it once
n_trial = 2;
R = CreateRecorder(O,60000); % The 2nd argument is Recorder.Capacity 
T_end = R.Capacity - 1; % simulation end time.  
%%
for i = 4:10 % Run 5 times
    Exp5;
    TransferSpikeRecord(R);
    BatchSTDPLearning(P_E);
    Load(R,1);
    Refresh(R);
end