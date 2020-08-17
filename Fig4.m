%% A localized initiation (Traveling wave is not supported at short range)
PanelFlag = 'A';
Exp4;
%% Prepare the figure
f = figure('Units','inches','Position',[1 0 4.25 8]);
ax = subplot(7,1,1);
for i = 2:7
    ax(i) = subplot(7,1,i);
end
set(ax,'Color','none','Box','off');
title(ax(1),'LVFA onset')
ylabel(ax(1),'Space');
ylabel(ax(2),'LFP (a.u.)');
title(ax(3),'Rhythmic onset');
ylabel(ax(3),'Space');
ylabel(ax(4),'LFP (a.u.)');
for i = [1:4,6]
    ax(i).XAxis.Visible = 'off';
end
ylabel(ax(5),{'Wave speed','cm/sec'});
ylabel(ax(6),{'Patient','LFP (mV)'});
ylabel(ax(7),{'Wave speed','cm/sec'});
set(ax,'NextPlot','add');
set(ax([2 4]),'YTick',[]);
set(ax([1 3]),'YTick',[]);
set(ax([1 3]),'YLim',[0 1]);
%%
TransferSpikeRecord(R);
plot(R);
l = R.Graph.Raster.Children.Children;
l.Parent = ax(1);
delete(R.Graph.Raster);
delete(R.Graph.Figure);
l.XData = l.XData / 1000; % Change unit to sec
l.YData = l.YData / max(O.n);
ax(1).YLim = [0 1];
pa_t = stim_t/1000;
pa_t = repmat(pa_t,[2 1]);
pa = patch(ax(1),pa_t(:)',[stim_x fliplr(stim_x)],[1 0 0],'FaceAlpha',0.3,'EdgeAlpha',0);
% clear LFP;LFP = LFP(R);
%% Period of interest
t_clip = stim_t/1000 + [-1 1];
t_clip = [t_clip;t_clip(2)+[20 25]];
%% Plot LFP
elec_position = [mean(stim_x) 0.33];
N_elec = round(elec_position * max(O.n));
C(1,:) = [0.47 0.67 0.19];
C(2,:) = [0.49 0.18 0.56];
for i = 1:numel(N_elec)
    l_LFP(i) = plot(ax(2),-LFP(N_elec(i),:),'LineWidth',0.5,'Color',C(i,:));
    l_LFP(i).XData = l_LFP(i).XData/1000;
end
l_LFP(2).YData = l_LFP(2).YData + max(abs(l_LFP(2).YData))*1.1;

% Clip the area you want to see
AxesClip(ax(1:2),t_clip,1);
axis(ax(1:2),'tight');
set(ax(2),'XLim',ax(1).XLim);
% After clip, you can add reference lines
for i = 1:numel(elec_position)
    l_neuron(i) = plot(ax(1),ax(1).XLim,[elec_position(i) elec_position(i)],'--','Color',C(i,:), 'LineWidth',0.5);
end



%% A hypersynchrnous initiation (Traveling wave is supported at short range
PanelFlag = 'B';
Exp4;
%% Put the result into the figure
TransferSpikeRecord(R);
plot(R);
l = R.Graph.Raster.Children.Children;
l.Parent = ax(3);
delete(R.Graph.Raster);
delete(R.Graph.Figure);
l.XData = l.XData / 1000; % Change unit to sec
l.YData = l.YData / max(O.n);
ax(3).YLim = [0 1];
pa_t = stim_t/1000;
pa_t = repmat(pa_t,[2 1]);
pa = patch(ax(3),pa_t(:)',[stim_x fliplr(stim_x)],[1 0 0],'FaceAlpha',0.3,'EdgeAlpha',0);
clear LFP;
LFP = LFP(R);
%% Period of interest
t_clip = stim_t/1000 + [-1 3];
t_clip = [t_clip;t_clip(2)+[15 18]];
%% Calculate wave speed
elec_position = [mean(stim_x) 0.31];
speed_position = elec_position(2);
N_elec = round(elec_position * max(O.n));
C(1,:) = [0.47 0.67 0.19];
C(2,:) = [0.49 0.18 0.56];
for i = 1:numel(N_elec)
    l_LFP(i) = plot(ax(4),-LFP(N_elec(i),:),'LineWidth',0.5,'Color',C(i,:));
    l_LFP(i).XData = l_LFP(i).XData/1000;
end
l_LFP(2).YData = l_LFP(2).YData + max(abs(l_LFP(2).YData))*1.1;

%% Calculate speed
N_speed = round((speed_position + [-0.015 0.015])* max(O.n));
sel = (R.S(2,:) >= N_speed(1))  &  (R.S(2,:) <= N_speed(2));
S = R.S(:,sel);
S(1,:) = S(1,:)/1000;
sel = false(1,size(S,2));
for i = 1:size(t_clip,1)
    sel = sel | (t_clip(i,1) < S(1,:) & t_clip(i,2) > S(1,:));
end
S = S(:,sel);
Tvec = t_clip(1):0.005:t_clip(end);
R_speed = histcounts(S(1,:),Tvec);
[pk,Idx_pk] = findpeaks(R_speed,'MinPeakProminence',20,'MinPeakDistance',10);
dT = [-0.07 0.07];
clear pValue Speed
for i = 1:numel(Idx_pk)
    sel = S(1,:)-Tvec(Idx_pk(i));
    sel = sel > dT(1) & sel < dT(2);
    S_local = S(:,sel);
    lm = fitlm(S_local(2,:)',S_local(1,:)');
    pValue(i) = lm.coefTest;    
    Speed(i) = 1/lm.Coefficients.Estimate(2);
end
S_speed = scatter(ax(5),Tvec(Idx_pk),abs(Speed),-log(pValue),'filled');
S_speed.YData = S_speed.YData/max(abs(S_speed.YData(pValue<0.001))); % normalize
S_speed.CData = (Speed'<0)*[1 0 1] + (Speed'>0)*[0.5 1 0.5];
% Clip the area you want to see
AxesClip(ax(3:5),t_clip,1);
axis(ax(1:5),'tight');
set(ax(4:5),'XLim',ax(3).XLim);
% After clip, you can add reference lines
% l_speed = plot(ax(3),ax(3).XLim,[speed_position speed_position],'--','Color',[0.85 0.33 0.10], 'LineWidth',0.5);
for i = 1:numel(elec_position)
    l_neuron(i) = plot(ax(3),ax(3).XLim,[elec_position(i) elec_position(i)],'--','Color',C(i,:), 'LineWidth',0.5);
end
l_speed_zero = plot(ax(5),ax(5).XLim,[0 0],'-','Color',[0.7 0.7 0.7], 'LineWidth',0.5);
ax(5).XTick = [];
ax(5).YLim = [0 1];
ax(5).YTick = 0;

%% Put in patient data 
m = matfile('D:\MATLAB\Eilepsy Codes\Method paper functions\method paper data\C5.mat');
lfp = m.LFP;
lfp_avg = mean(lfp);
[b,a]=butter(4,25/500); % Prevent ploting aliasing
lfp_avg = filtfilt(b,a,lfp_avg);
lfp_plot = plot(ax(6),(1:numel(lfp_avg))/1000,lfp_avg/1000,'Color',[0.3 0.75 0.93]);% change unit to mV
% ax(7)
m = matfile('D:\MATLAB\Eilepsy Codes\Method paper functions\RevisionRec.mat');
Rec =  m.RecMUA;
Rec = Rec(end,1);
Rec.S = sqrt(sum(Rec.V.^2,2));
Sc = scatter(ax(7),Rec.T(1:1049),Rec.S(1:1049),10,'filled');
Sc.SizeData = -log(Rec.p(1:1049));
CMap = hsv(360);
Color_idx = floor(rad2deg(Rec.direction+pi))+1;
Color_idx = Color_idx(1:1049);
Sc.CData =  CMap(Color_idx,:);
%%
t_clip = [318 338;378 388];
AxesClip(ax(6:7),t_clip,3);
%%
ax(6).Clipping = 'off';
axis(ax(6:7),'tight');
ax(6).YLim = [-1 1];
ax(6).YTick = 0;
ax(7).YLim = [0 70];
ax(7).YTick = [0 70];
ax(7).XTick = [];
ax(7).XLim = ax(6).XLim;
%% Add scale bar
ax_num = [1,3,6];
for i = 1:numel(ax_num)
    l_bar(i) = plot(ax(ax_num(i)),[1 2],[0.7 0.7],'LineWidth',2,'Color',[0.15 0.15 0.15]);
end
l_bar(3).XData = [2 3];
%% Add legend
colormap(ax(7),CMap);
caxis(ax(7),[-pi pi]);
cbar = colorbar('peer',ax(7),'North');
cbar.Limits = [-pi pi];
cbar.Ticks = [-pi pi];
cbar.TickLabels = {'-\pi','\pi'};
cbar.Box = 'off';
xlabel(cbar,'direction');
%% More legend
Sc_pValue = scatter(ax(5),[0.5 0.5 0.5],[0.5,0.7,0.9],[3,7,12],[0.15 0.15 0.15],'filled');
%%
Char = {'p=.05','p=10^-^3','p=10^-^5'};
for i = 1:3
    Tx_pValue(i) = text(ax(5),0.8,0.35 + 0.2*i,Char{i},'Color',[0.25 0.25 0.25]);
end
%%
Tx_dir = text(ax(5),5,0.4,'Forward','Color',[0.5 1 0.5]);
Tx_dir(2) = text(ax(5),9,0.4,'Backward','Color',[1 0 1]);

%%
Char = 'ABC';
for i = 1:numel(Char)
    annotation('TextBox','String',Char(i),'FontSize',14,'FontWeight','bold','EdgeColor','none')
end
%%
for i = 2:7
    ax(i).Position([1,3:4]) = ax(1).Position([1,3:4]);
end