% Fig 3 shows factors that affect pre-termination phase
% The experiment has been done in Fig2, so all you need is just to load it.

%% Fig 3A - zoom in
clear;close all;
F = openfig('D:\Paper Drafts\Model paper\Figures\2.fig');
rate = F.Children(end).Children(end).CData;
%%
f = figure('Units','inches','Position',[1 1 8 4]);
ax = axes(f,'Position',[0.07 0.7 0.9 0.25],'NextPlot','add');
T_pre_term = [67.2 77.2];
Im = imagesc(ax,T_pre_term,[0 1],rate(:,(T_pre_term(1)*1000):(T_pre_term(2)*1000)));
yo = ylabel(ax,'Space');
set(ax,'XLim',T_pre_term,'YLim',[0 1],'YDir','normal','XTIck',[],'YTick',[0 1]);
caxis(ax,[0 0.6]);
H = hot(360);
colormap(ax,H);
title(ax,'Pre-termination phase');
L_scale = plot([74.5 75.5],[0.85 0.85],'LineWidth',2,'Color',[0 1 1]);
Txt = text(75.8,0.85,'1 second','Color',[0 1 1]);
record_x = 0.3;
neuron_idx = round(500*record_x); 
L_stim = plot(T_pre_term,[record_x record_x],'LineStyle',':','Color',[0.4700    0.6700    0.1900],'LineWidth',1);

% Panel A - Lower subpanel Interburst interval
[Pk,Tb] = findpeaks(rate(neuron_idx,:),'MinPeakHeight',0.01);
Sel = T_pre_term(2) > Tb & Tb > T_pre_term(1);
Pk = Pk(Sel);
Tb = Tb(Sel);
ax_sub = axes(f,'Position',ax(1).Position,'NextPlot','add','Color','none','XTick',[],'YTick',[],'Clipping','off','YLim',[0 0.5],'XLim',T_pre_term);
ax_sub.Position(2) = 0.59;
ax_sub.Position(4) = 0.075;
ax_sub.YAxisLocation = 'right';
yo(2) = ylabel(ax_sub,'rate');
yo(2).Position(1) = T_pre_term(1) - 0.55;
rate_sub = rate(neuron_idx,(T_pre_term(1)*1000):(T_pre_term(2)*1000));
plot(ax_sub, T_pre_term(1):0.001:T_pre_term(2), rate_sub, 'Color', L_stim.Color,'LineWidth',1);
Tx_tick = text(ax_sub,T_pre_term(2)-0.8,0.5,'0.5 f_m_a_x','HorizontalAlignment','left');

% Panel B-E 
ax(2) = axes(f,'Position',[0.07 0.08 0.2 0.4]);
ax(3) = axes(f,'Position',[0.27+1/30 0.08 0.2 0.4]);
ax(4) = axes(f,'Position',[0.47+2/30 0.08 0.2 0.4]);
ax(5) = axes(f,'Position',[0.67+3/30 0.08 0.2 0.4]);
set(ax(2:5),'NextPlot','add','Color','none','YTick',[])
title(ax(2),'Firing rate');
title(ax(3),'Interburst interval');
title(ax(4),'Speed');
title(ax(5),'Burst width (RMSE)');
ylabel(ax(2),'spikes/sec');
ylabel(ax(3),'second');
ylabel(ax(4),'cm/sec');
ylabel(ax(5),'ms');

%% C5 Patient Data - this part requires you to connect the patient data hard drive
Ev = EventData.load('E:\data\Patient Data\C5\SeizureMUADataC503.mat');
map = [     -1    96    93    92    90    88    85    83    81    -1
            95    63    94    91    89    87    86    84    80    79
            32    61    59    57    55    53    49    82    78    77
            30    64    60    58    56    51    47    45    76    75
            28    31    62    52    46    44    43    41    74    73
            26    29    21    54    50    42    40    39    72    71
            24    27    25    19    15    48    38    37    70    69
            22    20    23    13    17     5    36    35    68    67
            18    16    12    11     9     7    34    33    66    65
            -1    14    10     8     6     4     3     1     2    -1];
ch = [     1     2     3     4     6     7     8     9    10    11    13    16    17    20    22    23    24 ...
    25    27    28    29    30    31    32    33    35    36    37    39    41    45    49    53    54 ...
    55    57    59    60    61    62    63    64    65    66    67    70    71    73    75    77    79 ...
    80    83    84    86    88    89    90    91    93    94    96]; % This is based on Cathy & Garret's UtahKey file
map(~ismember(map,ch)) = -1;
[~,Idx] = ismember(ch,map);
[p1,p2] = ind2sub(size(map),Idx);
%%
Tvec = 190:0.001:220; % The last 30 seconds
TS = Ev2TS(Ev,Tvec,'sigma',0.01); % We decrease kernel width to 10 ms to capture fast change
[~,pk] = findpeaks(sum(TS.Data(:,ch),2)/numel(ch),'minpeakheight',8,'minpeakdistance',100); % The same detection criteria as my previous paper
%%
for i = 1:numel(pk)
    Ev_local = select_event_time(Ev,Tvec(pk(i))+[-0.05 0.05],'include');
    T = Ev_local.timestamps(ch);
    P = cell(numel(ch),1);    
    for j=1:numel(ch)
        P{j} = repmat([p1(j),p2(j)], [numel(T{j}),1]);
    end
    LM = fitlm(cell2mat(P)*0.04,cell2mat(T)); % Change space unit to cm
    V(i,:) = pinv(LM.Coefficients.Estimate(2:3)); % cm/sec
    RMSE(i) = LM.RMSE;
    Speed(i) = norm(V(i,:));
    pValue(i) = LM.coefTest;
    A(i) = angle(V(i,1) + sqrt(-1)*V(i,2));
end
%%

% Peak firing rate
r = sum(TS.Data(:,ch),2)/numel(ch);
plot(ax(2),Tvec-Tvec(end),r,'Color',[0.6400    0.0800    0.1800]);
xlim(ax(2),[-15 0]);
ylim(ax(2),[0 200]);
ax(2).YTick = ax(2).YLim;
ax(2).XTick = ax(2).XLim;
ax(2).XTickLabel{1} = ['     ' num2str(ax(2).XLim(1)) 'sec'];
[ax(2).UserData.rho, ax(2).UserData.pValue] = corr(pk(pk>15000),r(pk(pk>15000)),'type','Spearman');

%% Interburst interval
sc = scatter(ax(3),Tvec(pk(2:end)) - Tvec(end),diff(Tvec(pk)),10,'filled');
xlim(ax(3),[-30 0]);
ylim(ax(3),[0 0.7]);
ax(3).YTick = ax(3).YLim;
ax(3).XTick = ax(3).XLim;
ax(3).XTickLabel{1} = ['     ' num2str(ax(3).XLim(1)) 'sec'];
[ax(3).UserData.rho, ax(3).UserData.pValue] = corr(sc.XData(:), sc.YData(:),'type','Spearman');

%% Speed
ax(4).ColorOrderIndex = 1;
sc = scatter(ax(4),Tvec(pk) - Tvec(end),Speed,1 + (pValue<0.05)*10,'filled');
xlim(ax(4),[-30 0]);
ylim(ax(4),[0 75]);
ax(4).YTick = ax(4).YLim;
ax(4).XTick = ax(4).XLim;
ax(4).XTickLabel{1} = ['     ' num2str(ax(4).XLim(1)) 'sec'];
[ax(4).UserData.rho, ax(4).UserData.pValue] = corr(sc.XData(sc.SizeData>1)', sc.YData(sc.SizeData>1)','type','Spearman');

%% RMSE (width)
ax(5).ColorOrderIndex = 1;
sc = scatter(ax(5),Tvec(pk) - Tvec(end),RMSE*1000,1 + (pValue<0.05)*10,'filled');
xlim(ax(5),[-30 0]);
ylim(ax(5),[0 25]);
ax(5).YTick = ax(5).YLim;
ax(5).XTick = ax(5).XLim;
ax(5).XTickLabel{1} = ['     ' num2str(ax(5).XLim(1)) 'sec'];
[ax(5).UserData.rho, ax(5).UserData.pValue] = corr(sc.XData(sc.SizeData>1)', sc.YData(sc.SizeData>1)','type','Spearman');

%% Save direction data in the figure
f.UserData.V = V;
f.UserData.A = A;
f.UserData.pValue = pValue;
f.UserData.T = Tvec(pk);
%%
Char = 'ABCDE';
for i = 1:numel(Char)
    annotation('TextBox','String',Char(i),'FontSize',14,'FontWeight','bold','EdgeColor','none')
end


%% Supplementary S5





%% Another way: directly load from old data ( This part is obsolete )
% load('D:\MATLAB\Eilepsy Codes\Method paper functions\method paper data\RevisionRec.mat','Record_MUA');
% Record_MUA(2) = [];
% sel = 950:1049; % The termination 100 burst
% Tvec = Record_MUA.T(sel);
F(2) = openfig('D:\Conference files\2016 Gordon Epilepsy & Synchronization\PosterFig4A.fig');
set(F(2).Children.Children,'Parent',ax(2));
F(3) = openfig('D:\Conference files\2016 Gordon Epilepsy & Synchronization\PosterFig4B.fig');
F(3).Children.Children.Parent = ax(3);
F(4) = openfig('D:\Conference files\2016 Gordon Epilepsy & Synchronization\PosterFig4C.fig');
F(4).Children.Children.Parent = ax(4);
F(5) = openfig('D:\Conference files\2016 Gordon Epilepsy & Synchronization\PosterFig4D.fig');
set(F(5).Children.Children,'Parent',ax(5));
delete(ax(5).Children(2)); % Remove the insignificant ones
close(F(2:5));
%%
XLIM = [85 100];
set(ax(2),'XLim',XLIM,'XTick',100,'XTickLabel',{[num2str(diff(XLIM)) ' sec    ']});
XLIM = [70 100];
set(ax(3:5),'XLim',XLIM,'XTick',100,'XTickLabel',{[num2str(diff(XLIM)) ' sec    ']});
ax(3).Children.MarkerSize = 10;
ax(4).Children.MarkerSize = 10;
ax(4).Children.YData = ax(4).Children.YData*1000;
ax(5).Children.MarkerSize = 10;
%%
ax(2).YLim = [0 200];ax(2).YTick = [0 200];
ax(3).YLim = [0 0.7];ax(3).YTick = [0 0.7];
ax(4).YLim = [0 30];ax(4).YTick = [0 30];
ax(5).YLim = [0 75];ax(5).YTick = [0 75];
ax(2).YLabel.Position(1) = 84;
ax(3).YLabel.Position(1) = 68;
ax(4).YLabel.Position(1) = 68;
ax(5).YLabel.Position(1) = 68;

%% Statistics
L = findobj(ax(2),'Type','Line');
L = L.YData(L.XData > 85 & L.XData < 100);
[pks,tpks] = findpeaks(L,'MinPeakProminence',10);
[Rho,pValue] = corr(tpks(:),pks(:),'type','Spearman')
disp(Rho)
disp(pValue)
%%
for i = 3:5
    L = findobj(ax(i),'Type','Line');
    X = L.XData(L.XData > 70 & L.XData < 100);
    Y = L.YData(L.XData > 70 & L.XData < 100);
    [Rho,pValue] = corr(X(:),Y(:),'type','Spearman')
end

%%
