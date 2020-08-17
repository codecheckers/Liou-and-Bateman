%% Figure 7 existence of spiral waves 
noise_sigma = 20; % Noisy environment
t_synchronization = 95000; % Timing of synchronization therapy
Exp7;close;

%% Prepare figure 
fig = figure('Units','inches','Position',[1 1 8 4.5]);
ax = axes('Position',[0.025,0.55,0.625 0.4],'Color','none','Box','off','NextPlot','add'); % all imagesc
set(ax,'YDir','reverse');
ax(1).XAxis.Visible = 'off';
ax(1).YAxis.Visible = 'off';
ax(2) = copyobj(ax(1),fig); % Spiral wave center panel
ax(2).Position = [0.7, ax(1).Position(2),0.2, ax(1).Position(4)]; 
ax(2).YDir = 'normal';
ax(3) = copyobj(ax(1),fig); % Successful stimulation
ax(3).Position = [ax(1).Position(1), 0.3, ax(1).Position(3), ax(1).Position(4)*10/32];
ax(4) = copyobj(ax(1),fig); % Failed stimulation 
ax(4).Position = [ax(1).Position(1), 0.1, ax(1).Position(3), ax(1).Position(4)*10/32];
ax(5) = copyobj(ax(1),fig); 
ax(5).Position = [ax(2).Position(1), ax(4).Position(2), ax(2).Position(3), ax(3).Position(2)+ax(3).Position(4)-ax(4).Position(2)];


%% Panel A (ax(1))
ax_idx = 1;
dx = fig.Position(3);
dy = fig.Position(4);
dt_vertical = 30000; % This determines how much time interval between the frames along vertical direction
ImIdx = 5000 + dt_vertical*(1:3);
% Vertical frames (coarse temporal evolution)
for i = numel(ImIdx):-1:1
    ImData = R.Var.V(:,ImIdx(i)) - R.Var.phi(:,ImIdx(i));
    ImData = reshape(ImData,O.n);
    ImData = O.param.f(ImData);
    ImData = ImData ./ O.param.f_max;
    ImData = flipud(ImData); % This is because the with imagesc, YDir of the axis is set to reverse by default
    ImData(~mask) = inf; % Applying the mask so that areas outside the boundar if white
    Im(i) = imagesc(ax(ax_idx),dx*([0 1]),dy*([0 1] + 1.1*(i-1)),ImData);
end

% Horizontal frames (fine temporal resolution)
IdxSel = [2 3];  % <----- this determines which coarse frames to zoom in temporally
dt_horizontal = 30; % This determines how much time interval between the frames along horizontal direction
for i = 1:numel(IdxSel)
    ImIdxSub = ImIdx(IdxSel(i))+dt_horizontal*(1:7); % 4 by 2 frames for each    
    for j = 1:numel(ImIdxSub) 
        ImData = R.Var.V(:,ImIdxSub(j)) - R.Var.phi(:,ImIdxSub(j));
        ImData = reshape(ImData,O.n);
        ImData = O.param.f(ImData);
        ImData = ImData ./ O.param.f_max;
        ImData = flipud(ImData); % This is because the with imagesc, YDir of the axis is set to reverse by default
        ImData(~mask) = inf; % Applying the mask so that areas outside the boundar if white
        ImSub(i,j) = imagesc(ax(ax_idx),dx*([0 1] + 1.1*j ),dy*([0 1] + 1.1*(IdxSel(i)-1)),ImData);    
    end
end    

axis(ax(1),'tight');
colormap(ax(1),hot(360));
caxis(ax(1),[0 0.6]);

cbar = colorbar('peer',ax(1),'South');
cbar.Position = [0.4 0.85 0.2 0.015];
cbar.AxisLocation = 'out';
cbar.Ticks = [0 0.6];
cbar.TickLabels = {'0','    0.6 f_m_a_x'};

%% Panel B (ax(2))
% Analysing center movement - see how many centers are there at the beginning
t0 = 75000; % <--------- Select some time when spiral waves have appeared so that we can trace the centers
v = R.Var.V(:,t0);
v(~mask(:)) = nan;
phi = R.Var.phi(:,t0);
phi(~mask(:)) = nan;
v_avg = median(v,'omitnan');
phi_avg = median(phi,'omitnan');
theta = angle(reshape((v-v_avg) + sqrt(-1)*(phi-phi_avg),O.n));
mask_reduced = mask;
mask_reduced(1:2,:) = 0;
mask_reduced(end-1:end,:) = 0;
mask_reduced(:,1:2) = 0;
mask_reduced(:,end-1:end) = 0;
mask_reduced = mask_reduced & ~conv2(~mask,ones(5),'same');
[Rmap,P0] = IdentifySpiralCenter(theta,'mask',mask_reduced); % So the P0 will be the initial centers to chase

% Trace each center to the end
t_search = t0:25:O.t; % <------ here it determines how often you are going to check where the centers are
I1 = cell(1,size(P0,1));
I2 = cell(1,size(P0,1));
for j = 1:size(P0,1)
    I1_seed = ceil(P0(j,1));
    I2_seed = ceil(P0(j,2));
    dI1 = -3:2; % Search range along dimension 1
    dI2 = -3:2; % Search range along dimension 2
    I1{j} = nan(1,numel(t_search));
    I2{j} = nan(1,numel(t_search));
    for i = 1:numel(t_search)
        v = R.Var.V(:,t_search(i));
        v(~mask(:)) = nan;
        phi = R.Var.phi(:,t_search(i));
        phi(~mask(:)) = nan;
        v_avg = median(v,'omitnan');
        phi_avg = median(phi,'omitnan');
        theta = angle(reshape((v-v_avg) + sqrt(-1)*(phi-phi_avg),O.n));
        [~,Posit] = IdentifySpiralCenter(theta,'Range1',I1_seed+dI1,'Range2',I2_seed+dI2,'mask',mask_reduced);
        if isempty(Posit)
            disp('The center has been annihilated.');break
        else
            I1{j}(i) = Posit(1);
            I2{j}(i) = Posit(2);
            I1_seed = ceil(Posit(1));
            I2_seed = ceil(Posit(2));
        end
    end
end

% Plotting the trajectories of the spiral wave centers
CMap = jet(numel(t_search));
for j = 1:size(P0,1)
    sel = I1{j} > 0;
    p1 = (I1{j}(sel)-1)/(O.n(1)-1);
    p1 = smooth(p1,41); % <----- here you need to decide how much smooth you want 
    p2 = (I2{j}(sel)-1)/(O.n(2)-1);
    p2 = smooth(p2,41); % <----- here you need to decide how much smooth you want 
    sc_path(j) = scatter(ax(2),p2,p1,1,CMap(1:numel(p1),:));hold on;
end % PS: The reason why smoothing is required is because the state space theta estimate is noisy due to I_noise

% Put on the outside circle
theta = 0:0.01:(2*pi);
l_circle = plot(ax(2),0.5 + 0.5*cos(theta),0.5 + 0.5*sin(theta),'Color',[0.15 0.15 0.15]);

% Time color bar
TLim = t_search/1000; % Time unit: sec
caxis(ax(2),[min(TLim),max(TLim)]);
colormap(ax(2),jet(numel(p1)));
cbar_time = colorbar('peer',ax(2),'east');
cbar_time.Position = [0.93 ax(2).Position(2) 0.015 ax(2).Position(4)];
cbar_time.Limits = [75 200];
cbar_time.Ticks = cbar_time.Limits;
cbar_time.TickLabels = cellfun(@(x) [x ' s'], cbar_time.TickLabels,'Un',0);
cbar_time.AxisLocation = 'out';
title(ax(2),'Spiral centers');
axis(ax(2),'equal');

%% Panel C & D (ax(3) & ax(4)): try synchronization therapy 
% First keep the original simulation data safe
% OriginalRecord = R.Var;

for idx_ax = 4:-1:3 % Index of axes, the reason for doing 3 later is for making videos
    % Load the time when the last coarse frame was shown
    T_load = ImIdx(end);
    Load(R,T_load); 
    O.Input.E = O.UserData.Input.E;
    O.Input.I = O.UserData.Input.I;
    O.Proj.In(1).Value = O.UserData.Proj.In(1).Value;
    O.Proj.In(2).Value = O.UserData.Proj.In(2).Value;
    O.Proj.In(3).Value = O.UserData.Proj.In(3).Value;

    % Build a shock external input (successful pulse therapy)
    O.Ext(end+1) = ExternalInput;
    O.Ext(end).Target = O;
    switch idx_ax
        case 3
            A_shock = 200; % Successful treatment
        case 4
            A_shock = 100; % Failed treatment 
    end
    T_shock = 30;
    O.Ext(end).Deterministic = @(x,t) A_shock * ones(O.n); % x: position, t: ms, unit: mV
    O.Ext(end).Tmax = O.t + T_shock; % 100 ms;

    for i = 1:T_shock+3000
        Update(O,1);
        WriteToRecorder(O);     
    end

    % Plot the results into ax(3)
    ImIdxSub = ImIdx(end) + dt_horizontal*(1:7);
    for j = 1:numel(ImIdxSub) 
        ImData = R.Var.V(:,ImIdxSub(j)) - R.Var.phi(:,ImIdxSub(j));
        ImData = reshape(ImData,O.n);
        ImData = O.param.f(ImData);
        ImData = ImData ./ O.param.f_max;
        ImData = flipud(ImData);
        ImData(~mask) = inf;
        imagesc(ax(idx_ax),...
                dx*([0 1] + 1.1*j ), ...
                dy*([0 1]), ...
                ImData);
    end
end
axis(ax(3:4),'tight');
caxis(ax(3),ax(1).CLim);
caxis(ax(4),ax(1).CLim);
colormap(ax(3),hot(360));
colormap(ax(4),hot(360));
ax(3).XLim = ax(1).XLim;
ax(4).XLim = ax(1).XLim;
% % Plot a pulse
% l_stim = plot(ax(3),1+[0 2 2 4 4 6],3-[0 0 2 2 0 0],'Color',[0.15 0.15 0.15]);
% Tx_stim = text(ax(3),0,0,[num2str(T_shock) 'ms']);
% Tx_stim(2) = text(ax(3),0,0,[num2str(A_shock) 'pA'],'Rotation',90);

%% Panel E - Try different shock amplitude/duration 
result_known = 1; % This is just a switch so that if you have the results already, just load it.
if result_known
    fig_supp = openfig('D:\Paper Drafts\Model paper\7Extra.fig');
    copyobj(fig_supp.Children(2).Children,ax(5));
else
    % Build a shock external input
    A_shock = 25:25:300;
    T_shock = 2:2:50;

    for i = numel(A_shock):-1:1
        for j = numel(T_shock):-1:1
            Load(R,t_synchronization); % This one only loads the main variables ... (The object-oriented code needs to be improved in the future)    
            O.Input.E = O.UserData.Input.E;
            O.Input.I = O.UserData.Input.I;
            O.Proj.In(1).Value = O.UserData.Proj.In(1).Value;
            O.Proj.In(2).Value = O.UserData.Proj.In(2).Value;
            O.Proj.In(3).Value = O.UserData.Proj.In(3).Value;

            O.Ext(end+1) = ExternalInput;
            O.Ext(end).Target = O;          
            O.Ext(end).Deterministic = @(x,t) A_shock(i) * ones(O.n); % x: position, t: ms, unit: mV
            O.Ext(end).Tmax = O.t + T_shock(j); % After this, the ExternalInput object will be deleted

            for iter = 1:T_shock(j)+3000
                Update(O,1);

                if iter == T_shock(j)+1
                    % Immediately after the shock
                    V_distribution{i,j} = O.V(mask);
                    phi_distribution{i,j} = O.phi(mask);
                    Cl_in_distribution{i,j} = O.Cl_in(mask);                
                    g_K_distribution{i,j} = O.g_K(mask);                
                end
            end

            % Mean firing rate 3 seconds after the shock therapy
            f_end = O.param.f(O.V-O.phi);
            f_distribution{i,j} = f_end(mask);
        end
    end
    imagesc(ax(5), T_shock, A_shock, cellfun(@max,f_distribution)/0.2); % Because max rate is 0.2 kHz
end
axis(ax(5),'tight');
ax(5).YDir = "normal";
ax(5).YTick = 300;
ax(5).YAxis.Visible = 'on';
ax(5).XTick = 50;
ax(5).XTickLabel = {'   50 ms'};
ax(5).XAxis.Visible = 'on';

colormap(ax(5),jet(120));
cbar_pulse = colorbar('peer',ax(5),'east');
cbar_pulse.Position = [cbar_time.Position(1), ax(5).Position(2), cbar_time.Position(3), ax(5).Position(4)];
cbar_pulse.Limits = [0 0.55];
cbar_pulse.Ticks = cbar_pulse.Limits;
cbar_pulse.AxisLocation = 'out';

title(ax(5),'Syn. pulse effectiveness')

%%
export_fig 7 -png -transparent -r600 -c[0,0,0,0];
%% Create a video for demonstration
MakeVideo = true;% Just a switch

if MakeVideo 
    VIDEO = VideoWriter('Figure7Video','MPEG-4');
    VIDEO.Quality  = 55;
    VIDEO.FrameRate = 100;
    open(VIDEO)
    f = figure('Color',[1 1 1]);f.Units = 'inches';f.Position = [1 1 1.9 1.5];
    ax_video = axes('Position',[0 0.05 1 0.9]);
    ImData = O.param.f(reshape(R.Var.V(:,1) - R.Var.phi(:,1),O.n(1),O.n(2)))./O.param.f_max;
    ImData(ImData == 0) = inf;
    Im = imagesc(ax_video,ImData);
    ax_video.XAxis.Visible = 'off';
    ax_video.YAxis.Visible = 'off';
    ax_video.Units = 'inches';
    ax_video.Position(3) = ax_video.Position(4);
    ax_video.Units = 'normalized';
    ax_video.YDir = 'normal';
    tx_video = text(ax_video,53,4,{'0.00','sec'},'HorizontalAlignment','center');
    caxis([0 0.6])
    H = hot(360);
    colormap(ax_video,H);
    cbar_video = colorbar('peer',ax_video,'east');
    cbar_video.Position = [0.84 0.05 0.05 0.9];
    cbar_video.Ticks = [0 0.6];
    cbar_video.TickLabels{2} = '.6';
    yo = ylabel(cbar_video,'f_m_a_x');
    yo.Position(1) = 1;
    set(findall(f,'Units','normalized'),'Units','inches');
    f.Position(4) = 1.9;   
    ax_video_title = title(ax_video,{'', ...
                                     ''});
    ax_video_title.Position(1) = 30;
    
    for i = 10:10:ImIdx(end) % Here is the pre-stimulation part
        if i == 2000
            ax_video_title.String = 'Seizure onset';
        elseif i == 6000
            ax_video_title.String = '';
        elseif i == 10000
            ax_video_title.String = {'Tonic-to-clonic'; ...
                                     '  transition  '};
        elseif i == 15000
            ax_video_title.String = '';            
        elseif i == 40000
            ax_video_title.String = {'Spiral wave', ...
                                     ' emergence '};
        elseif i == 45000
            ax_video_title.String = '';
        elseif i == 60000
            ax_video_title.String = {'Spiral wave', ...
                                     'interactions'};
        elseif i == 70000
            ax_video_title.String = '';       
        elseif i == 77000
            ax_video_title.String = {'Status', ...
                                     'epilepticus'};
        elseif i == 92000
            ax_video_title.String = '';                
        end
        ImData = O.param.f(reshape(R.Var.V(:,i) - R.Var.phi(:,i),O.n(1),O.n(2)))/0.2;
        ImData(ImData == 0) = inf;
        Im.CData = ImData;
        drawnow;
        im = getframe(gcf);
        tx_video.String = {num2str(i/1000,'%.2f'),'sec'};
        writeVideo(VIDEO,im);
    end
    
    ax_video_title.String = {'Synchronizing pulse', ...
                             '   30 ms, 200 pA   '};drawnow;
    im = getframe(gcf);                         
    
    for i = 1:500 % Give 5 second pause in between for the reader to react
        writeVideo(VIDEO,im);        
    end
    
    for i = ImIdx(end) + (1:1000) % Here is the post-stimulation part
        ImData = O.param.f(reshape(R.Var.V(:,i) - R.Var.phi(:,i),O.n(1),O.n(2)))/0.2;
        ImData(ImData == 0) = inf;
        Im.CData = ImData;
        drawnow;
        im = getframe(gcf);
        tx_video.String = {num2str(i/1000,'%.2f'),'sec'};
        if i == ImIdx(end) + 300
            ax_video_title.String = 'Termination';              
        end
        writeVideo(VIDEO,im);
    end
    
    close(VIDEO);
end
