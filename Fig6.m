%% Prepare the figure
fig = figure('Units','inches','Position',[1 1 8 5]);
ax = axes('Position',[0.035,0.1,0.25*10/11 0.8],'Color','none','Box','off','NextPlot','add'); % all imagesc
ax(2) = copyobj(ax(1),fig); % Noisy version
ax(3) = copyobj(ax(1),fig); % Peak samples
ax(4) = copyobj(ax(1),fig); % Peak meet 
ax(5) = copyobj(ax(1),fig); % Direction
ax(2).Position(1) = 0.315;
ax(3).Position = [0.62 0.8 0.35 0.1];
ax(3).Position(4) = ax(3).Position(3) * 1/6 * 8 / 5;
ax(3).Position(2) = 0.9 - ax(3).Position(4);
ax(4).Position = [0.62 0.5 0.35 0.3];
ax(5).Position = [0.62 0.1  0.35 0.27];

set(ax(1:3),'YDir','reverse');
ylim(ax(5),[-pi pi]);
ax(5).YTick = [-pi 0 pi];
ax(5).YTickLabel = {'-\pi' '0' '\pi'};

ax(1).XAxis.Visible = 'off';
ax(2).XAxis.Visible = 'off';
ax(1).YAxis.Visible = 'off';
ax(2).YAxis.Visible = 'off';

%% Simulation & Figure 6A, 6B
ax_idx = 1;
for noise_sigma = [0 20]
    % Do the experiment
    Exp6;close;
    % Plot its results to ax(1) & ax(2) respectively, depending on ax_idx
    dx = fig.Position(3);
    dy = fig.Position(4);
    % Horizontal: coarse time scale evolution (the phasic evolution)
    ImIdx = 8500 + 10000*(0:4) + 305; % <------ this determines which frames to plot
    for i = numel(ImIdx):-1:1
        ImData = R.Var.V(:,ImIdx(i)) - R.Var.phi(:,ImIdx(i));
        ImData = reshape(ImData,O.n);
        ImData = O.param.f(ImData);
        ImData = ImData ./ O.param.f_max;
        ImData = flipud(ImData); % This is because default axes where imagesc is used 'YDir' is 'reverse'
        ImData(~mask) = inf; % Apply boundary mask, so that the image outside are all 'white'
        Im(i) = imagesc(ax(ax_idx),dx*([0 1] + 1.1*(i-1)),dy*([0 1]),ImData);
    end
    % Vertical: fine time scale evolution (back traveling waves in 2D)
    IdxSel = [1 3 5]; % <----- this determines which coarse frames to zoom in temporally
    for i = 1:numel(IdxSel)
        dt_frame = 20;
        ImIdxSub = ImIdx(IdxSel(i))+dt_frame*(1:10); % This determines how 'fine' the vertical frames are temporally
        for j = 1:numel(ImIdxSub) 
            ImData = R.Var.V(:,ImIdxSub(j)) - R.Var.phi(:,ImIdxSub(j));
            ImData = reshape(ImData,O.n);
            ImData = O.param.f(ImData);
            ImData = ImData ./ O.param.f_max;
            ImData = flipud(ImData); % This is because default axes where imagesc is used 'YDir' is 'reverse'
            ImData(~mask) = inf; % Apply boundary mask, so that the image outside are all 'white'
            ImSub(i,j) = imagesc(ax(ax_idx),dx*([0 1] + 1.1*(IdxSel(i)-1)),dy*([0 1] + 1.1*j),ImData);    
        end
    end
    caxis(ax(ax_idx),[0 0.6]);
    colormap(ax(ax_idx),hot(360));
    axis(ax(ax_idx),'tight');

    
    %% Here you can also save the results into a video
    make_video = true; % <----------------- Here you decide if you want to make video, too.
    if make_video 
        dN = 10; % Frame rate: 10 ms/frame
        VIDEO = VideoWriter(['Figure6Video' num2str(ax_idx)],'MPEG-4');
        VIDEO.Quality  = 55;
        VIDEO.FrameRate = 100;
        open(VIDEO)
        f_video = figure('Color',[1 1 1]);f_video.Units = 'inches';f_video.Position = [1 1 1.75 1.35];
        ax_video = axes('Position',[0 0.05 1 0.9]);
        rate = (R.Var.V - R.Var.phi);
        ImDataVideo = O.param.f(reshape(rate(:,1),O.n(1),O.n(2)))/ O.param.f_max(1);
        ImDataVideo(ImDataVideo == 0) = inf;
        Im = imagesc(ax_video,ImDataVideo);
        ax_video.XAxis.Visible = 'off';
        ax_video.YAxis.Visible = 'off';
        ax_video.Units = 'inches';
        ax_video.Position(3) = ax_video.Position(4);
        ax_video.Units = 'normalized';
        tx_video = text(ax_video,O.n(1)+3.5,O.n(2)-4.5,{'0.00','sec'},'HorizontalAlignment','center','FontSize',9);
        caxis([0 0.5])
        H = hot(360);
        colormap(ax_video,H);
        cbar_video = colorbar('peer',ax_video,'east');
        cbar_video.Position = [0.84 0.05 0.05 0.9];
        cbar_video.Ticks = [0 0.5];
        cbar_video.TickLabels{2} = '.5';
        yo = ylabel(cbar_video,'f_m_a_x');
        yo.Position(1) = 1;

        for i = 1:7000
            ImDataVideo = O.param.f(reshape(rate(:,1+dN*i),O.n(1),O.n(2)))/O.param.f_max(1);
            ImDataVideo(ImDataVideo == 0) = inf;
            Im.CData = ImDataVideo;
            drawnow;
            im = getframe(gcf);
            tx_video.String = {num2str(i/100,'%.2f'),'sec'};
            writeVideo(VIDEO,im);
        end
        close(VIDEO);         
    end
    %%
    ax_idx = ax_idx + 1;    
end
title(ax(1),'Noise free');
title(ax(2),['White noise, \sigma=' num2str(noise_sigma) 'pA']);
% the associated colorbar
cbar = colorbar('peer',ax(1),'South');
cbar.Position = [0.05 0.065 0.2 0.01];
cbar.AxisLocation = 'out';
cbar.Ticks = [0 0.6];
cbar.TickLabels = {'0','    0.6 f_m_a_x'};

%% Analyse traveling weaves & peak points
% Wave calculation settings 
% find wave recording electrodes
[Py,Px] = meshgrid(1:O.n(1),1:O.n(2));
P = [Px(:),Py(:)]; % Neuronal 2d indices
record_x = 0.5; % Normalized spatial unit
record_y = 0.5; % Normalized spatial unit
record_r = 0.05; % Normalized spatial unit
record_sel = (sqrt(((P(:,1)-1)/(O.n(2)-1)-record_x).^2 + ((P(:,2)-1)/(O.n(1)-1)-record_y).^2)) < record_r;
P = P(record_sel,:);

% Confluence point chasing settings
n_center = ceil(O.n/2); 
idx_center = sub2ind(O.n,n_center(1),n_center(2)); % We chase from the center of the field
dT_max = 10; % maximally it can have 10 ms delay in terms of peak timings when 'chasing confluence points' 
P_confluence = [];

% Find bursts for eahc channel (oscillatory voltage peaks)
t_start = 8000; % Since this time you start to look for bursts
v = R.Var.V(:,t_start+1:end) - R.Var.phi(:,t_start+1:end);
v_pk = cell(size(v,1),1);
t_pk = cell(size(v,1),1);
for i = size(v,1):-1:1
    [v_pk{i},t_pk{i}] = findpeaks(v(i,:),'minpeakdistance',100,'minpeakprominence',5);
end
t_pk = cellfun(@(x) x + t_start,t_pk,'un',false);

% Loop to chase confluence point & estimate wave direction
for i = 1:numel(t_pk{idx_center})
    T = [];
    Tref = t_pk{idx_center}(i);
    % construct the delta_T map with respective to the center 
    for j = numel(t_pk):-1:1 % j : index of channel
        if isempty(t_pk{j})
            T(j) = nan;
        else
            dT = t_pk{j}-Tref;
            [~,idx] = min(abs(dT));
            T(j) = dT(idx);
        end
    end
    % Calculate wave statistics
    T_for_wave = T(record_sel);
    T_for_wave(abs(T_for_wave) > 100) = nan; % If the peak is too temporally distant, discard it
    if sum(~isnan(T_for_wave))<=3
        V(i,:) = [nan nan];
        Speed(i) = nan;
        A(i) = nan;
        F_test(i) = nan;
    else
        LM = fitlm(P,T_for_wave);
        V(i,:) = pinv(LM.Coefficients.Estimate(2:3)); % Velocity
        Speed(i) = norm(V(i,:)); % Speed
        A(i) = angle(V(i,2) + sqrt(-1)*V(i,1)); % Wave traveling direction
        F_test(i) = LM.coefTest; % Linear regression F-test pValue
    end
    % reshape it
    T = reshape(T,O.n);
    % Chase peak
    P_confluence(i,:) = ChasePeak(T,n_center);
end

%% Panel C, show a few sample episodes' confluence points
sample_idx = 90:95; % <----- here determines which ones to show
sc_peak = scatter(ax(3),(1:numel(sample_idx))'*50+P_confluence(sample_idx,2),P_confluence(sample_idx,1),8,'filled');hold on
sc_peak.CData = [0.9300    0.6900    0.1300];
sc_center = scatter(ax(3),(1:numel(sample_idx))*50 + 25, (1:numel(sample_idx))*0+25,'+');
sc_center.CData = [0.15 0.15 0.15];
theta = 0:0.01:(2*pi);
for i = 1:numel(sample_idx)
    l_circle(i) = plot(ax(3),25+50*i+25*cos(theta),25 + 25*sin(theta),'Color',[0.15 0.15 0.15]);
end
ax(3).XAxis.Visible = 'off';
ax(3).YAxis.Visible = 'off';
title(ax(3),'Confluence spot');

%% Panel D, show confluence point evolution
ax(4).ColorOrderIndex = 5;
l_peak = plot(ax(4),(P_confluence-1)./(O.n-1) - 0.5,'.');
axis(ax(4),'tight');
ax(4).YLim = [-0.5 0.5];
ax(4).XAxis.Visible = 'off';
ax(4).YTick = [-0.5 0.5];
phase_mark = 120; % <---- here determines when wavefront is annihilated
lm = plot(ax(4),[phase_mark phase_mark],ax(4).YLim,'--','Color',[0.5 0.5 0.5]);
lm2 = plot(ax(4),sample_idx(1),0.35,'v','Color',sc_peak.CData);
lm2.MarkerFaceColor = sc_peak.CData;
le = legend(l_peak,'y-coordinate','x-coordinate');
le.Color = 'none';
le.EdgeColor = 'none';
Tx = text(ax(4),phase_mark/2,0.4,{'Ictal','clonic'},'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
Tx(2) = text(ax(4),mean([phase_mark ax(4).XLim(2)]),0.4,{'Pre','termination'},'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
yo = ylabel(ax(4),'Space');
yo.Position = [-6    0   -1];

%% Panel E, traveling wave directions
sca_A = scatter(ax(5),1:numel(A),A,-log(F_test),'filled');
lm3 = plot(ax(5),[phase_mark phase_mark],ax(5).YLim,'--','Color',[0.5 0.5 0.5]);
ax(5).XLim = ax(4).XLim;
title(ax(5),'Wave direction at the center');
Sc_pValue = scatter(ax(5),[40 40 40],[-2,-1,0],[3,7,12],[0.15 0.15 0.15],'filled');
Char = {'p=.05','p=10^-^3','p=10^-^5'};
for i = 1:3
    Tx_pValue(i) = text(ax(5),50,-3 + i+0.2,Char{i},'Color',[0.25 0.25 0.25],'VerticalAlignment','middle');
end
Ta = annotation('textarrow',[0.38 0.41],[0.065 0.065],'String','10 sec/frame');
Ta(2) = annotation('arrow',[0.54 0.54],[0.09 0.03]);
Ta_string = annotation('textbox',[0.426,0.04,0.12,0.06],'String','20 ms/frame','EdgeColor','none');
xlabel(ax(5),'Discharge index');

%% Calculate statistics
disp('Ydirection')
[rho,pValue] = corr(P_confluence(1:phase_mark-1,1),P_confluence(2:phase_mark,1));
fig.UserData.rho.Y = rho;
fig.UserData.p.Y = pValue;
disp('Xdirection')
[rho,pValue] = corr(P_confluence(1:phase_mark-1,2),P_confluence(2:phase_mark,2));
fig.UserData.rho.X = rho;
fig.UserData.p.X = pValue;

Asub = A;
Asub(F_test > 0.001) = nan;
Asub1 = Asub(1:phase_mark);
Asub1(isnan(Asub1)) = [];
[rho,pValue] = circ_corrcc(Asub1(1:end-1), Asub1(2:end));
fig.UserData.rho.dir.clonic = rho;
fig.UserData.p.dir.clonic = pValue;

Asub1 = Asub(phase_mark+1:end);
Asub1(isnan(Asub1)) = [];
[rho,pValue] = circ_corrcc(Asub1(1:end-1), Asub1(2:end));
fig.UserData.rho.dir.preterm = rho;
fig.UserData.p.dir.preterm = pValue;

%% Panel labels
Char = 'ABCD';
for i = 1:numel(Char)
    An(i) = annotation('TextBox',[0.1 0.9 0.1 0.1],'Unit','normalized','String',Char(i),'FontSize',14,'FontWeight','bold','EdgeColor','none');
end

