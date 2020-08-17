%% Fig 1 - generating video and show the similarity between data & model 
%Exp1;

%% Prepare figure
fig = figure('Units','inches','Position',[1 0.5 4.5 7.5]);

% Slow time scale images
ax = axes('Position',[0.07 0.08 0.15 0.6],'NextPlot','add','Color','none','Box','off', ...
          'XTick',[],'YTick',[],'YDir','reverse'); % Slow time scale, data
n_img = 6; % How many images per column
image_ratio = (7+1)/9;
ax.Units = 'inches';
ax.Position(4) = ax.Position(3)*n_img*image_ratio;
ax.Units = 'normalized';
      
ax(2) = copyobj(ax,fig); % Slow time scale, simulation results
ax(2).Position(1) = 0.255;
ax(2).YAxisLocation = 'right';

% Fast time scale imagesc
ax(3) = copyobj(ax(1),fig);
ax(3).Position(1) = 0.6;
ax(4) = copyobj(ax(1),fig);
ax(4).Position(1) = 0.7855;
ax(4).YAxisLocation = 'right';

ax_model = copyobj(ax(1),fig);
ax_model.Position(2) = ax(1).Position(2) + ax(1).Position(4) + 0.17;
ax_model.Position(3) = ax(4).Position(1) + ax(4).Position(3) - ax_model.Position(1);
ax_model.Position(4) = 0.18;

%% Patient data part (ax(1) & ax(3))
% Slow time scale: ax(1)
% m = matfile('Data/elife-50927-data1-v1.mat'); 
% map = m.map; % data map
% rate = m.R100; % Slow time scale data
% rate = InterchangeMatrixArray(rate.Data','map',map);
% n_start = 200;
% for i = 1:n_img
%     imagesc(ax(1),[0 1],[0 7/9] + (i-1)*1,(rate(:,:,n_start+i*100)) ); % frame rate : 1/sec
% end
% ylabel(ax(1),'1 sec/frame');
% CLim = [0 300];
% caxis(CLim);
H = hot(360);
% colormap(ax(1),H);
% cbar = colorbar('peer',ax(1),'south');
% cbar(1).Position = ax(1).Position;
% cbar(1).Position(2) = cbar(1).Position(2) - 0.02;
% cbar(1).Position(4) = 0.01;
% cbar(1).Ticks = CLim;
% cbar(1).TickLabels{1} = ['        ' num2str(CLim(1)) ' spks/sec'];
% cbar(1).TickLabels{2} = [num2str(CLim(2)) '  '];
% xlabel(cbar(1),'(100ms kernel)');
% Ta = annotation('TextArrow',repmat(ax(1).Position(1)-0.01,[2 1]), ...
%                             [ax(1).Position(2)+ax(1).Position(4), ax(1).Position(2)], ...
%                             'LineWidth',1);
% 
% % Fast time scale - ax(3)
% n_start = 9198;
% rate = m.R10; % Slow time scale data
% rate = InterchangeMatrixArray(rate.Data','map',map);
% for i = 1:n_img
%     imagesc(ax(3),[0 1],[0 7/9] + (i-1),(rate(:,:,n_start+i*10))); % every 10 ms 
% end
% ylabel(ax(3),'10 ms/frame');
% CLim = [0 300];
% caxis(ax(3),CLim);
% colormap(ax(3),H);
% cbar(3) = colorbar('peer',ax(3),'south');
% cbar(3).Position = ax(3).Position;
% cbar(3).Position(2) = cbar(3).Position(2) - 0.02;
% cbar(3).Position(4) = 0.01;
% cbar(3).Ticks = CLim;
% cbar(3).TickLabels{1} = ['        ' num2str(CLim(1)) ' spks/sec'];
% cbar(3).TickLabels{2} = [num2str(CLim(2)) '  '];
% 
% xlabel(cbar(3),'(10ms kernel)');
% Ta(3) = annotation('TextArrow',repmat(ax(3).Position(1)-0.01,[2 1]), ...
%                               [ax(1).Position(2)+ax(1).Position(4), ax(1).Position(2)], ...
%                               'LineWidth',1);



%% Simulation slow time scale data
n_start = 1100;
window_x = 8 + (1:27);
window_y = 33 + (1:21);
for i = 1:n_img
    rate = 0;
    for j = 1:20 % 200 ms avg
        ij = n_start + i*200 + j;
        v = R.Var.V(:,ij) - R.Var.phi(:,ij);
        v = reshape(v,O.n(1),O.n(2),[]);        
        rate = rate + O.param.f(v);
    end
    rate = (rate / j) / mean(O.param.f_max(:));
    imagesc(ax(2),[0 1],[0 7/9] + (i-1)*0.9,rate(window_y,window_x)); % every 10 ms 
end

ylabel(ax(2),'2 sec/frame');
CLim = [0 0.22];
caxis(ax(2),CLim);
colormap(ax(2),H);
cbar(2) = colorbar('peer',ax(2),'south');
cbar(2).Position = ax(2).Position;
cbar(2).Position(2) = cbar(2).Position(2) - 0.02;
cbar(2).Position(4) = 0.01;
cbar(2).Ticks = CLim;
cbar(2).TickLabels{2} = ['    ' num2str(CLim(2)) ' f_m_a_x'];
xlabel(cbar(2),'(200ms avg)');
Ta(2) = annotation('TextArrow',repmat(ax(2).Position(1)+ax(2).Position(3)+0.01,[2 1]), ...
                            [ax(2).Position(2)+ax(2).Position(4), ax(2).Position(2)], ...
                            'LineWidth',1);
                        
% Fast time scale
n_start = 2995;
for i = 1:n_img
    v = R.Var.V(:,n_start + 2*i) - R.Var.phi(:, n_start + 2*i); % every 20 ms = 1 frame
    v = reshape(v,O.n(1),O.n(2),[]);
    rate = (O.param.f(v)/mean(O.param.f_max(:)));
    imagesc(ax(4),[0 1],[0 7/9] + (i-1)*0.9,rate(window_y,window_x)); % every 10 ms 
end

ylabel(ax(4),'20 ms/frame');
CLim = [0 0.5];
caxis(ax(4),CLim);
colormap(ax(4),H);
cbar(4) = colorbar('peer',ax(4),'south');
cbar(4).Position = ax(4).Position;
cbar(4).Position(2) = cbar(4).Position(2) - 0.02;
cbar(4).Position(4) = 0.01;
cbar(4).Ticks = CLim;
cbar(4).TickLabels{2} = ['    ' num2str(CLim(2)) ' f_m_a_x'];
xlabel(cbar(4),'(Instantaneous)');
Ta(4) = annotation('TextArrow',repmat(ax(4).Position(1)+ax(4).Position(3)+0.01,[2 1]), ...
                            [ax(4).Position(2)+ax(4).Position(4), ax(4).Position(2)], ...
                            'LineWidth',1);


%% Aftermath
axis(ax,'tight');
Ot = title(ax(1),'Patient');Ot.Position(2) = -1.35;
Ot(3) = title(ax(3),'Patient');Ot(3).Position(2) = -1.35;
Ot(2) = title(ax(2),'Model');Ot(2).Position(2) = -1.17;
Ot(4) = title(ax(4),'Model');Ot(4).Position(2) = -1.17;

for i = 1:numel(ax)
    ax(i).XAxis.Color = [1 1 1];
    ax(i).YAxis.Color = [1 1 1];
    ax(i).YAxis.Label.Color = [0.15 0.15 0.15];
end

%% load the schematic
% p = imread('D:\Paper Drafts\Model paper\Figures\1A.png', 'png', 'BackgroundColor',[1 1 1]);
% ImSchematic = image(ax_model, p);
% axis(ax_model,'image');
% axis(ax_model,'manual');
% ax_model.XTick = [];
% ax_model.YColor = [1 1 1];
% ax_model.Position(3:4)=ax_model.Position(3:4)*1.4 ;
% TxModel = text(ax_model,0,0,'Excitatory','Color',[1 0 0],'FontWeight','bold');
% TxModel(2) = text(ax_model,0,0,'Inhibitory','Color',[0 0 1],'FontWeight','bold');

%% Add schematic illustrations
ax_map = axes('Position',[ax(1).Position(1), ax(1).Position(2)+ax(1).Position(4)+0.005, ax(1).Position(3), 0.09], ...
              'Color','none','XTick',[],'YTick',[],'Box','off','NextPlot','add', ...
              'XColor',[1 1 1],'YColor',[1 1 1],'YDir','reverse');
ax_map.Units = 'inches';
ax_map.Position(4) = ax_map.Position(3);
ax_map.Units = 'normalized';
ax_map(2) = copyobj(ax_map(1),fig);
ax_map(2).Position(1) = ax(2).Position(1);
ax_map(3) = copyobj(ax_map(1),fig);
ax_map(3).Position(1) = ax(3).Position(1);
ax_map(4) = copyobj(ax_map(1),fig);
ax_map(4).Position(1) = ax(4).Position(1);

theta = 0:0.01:2*pi;
plot(ax_map(2),0.5 + 0.5*cos(theta), 0.5+0.5*sin(theta),'Color',[0.5 0.5 0.5],'LineWidth',1);
plot(ax_map(4),0.5 + 0.5*cos(theta), 0.5+0.5*sin(theta),'Color',[0.5 0.5 0.5],'LineWidth',1);

x_sq = [window_x(1),window_x(end)]/O.n(2);
y_sq = [window_y(1),window_y(end)]/O.n(1);
plot(ax_map(2),[x_sq,fliplr(x_sq),x_sq(1)],[y_sq(1),y_sq(1),y_sq(2),y_sq(2),y_sq(1)],'Color',[0.4700    0.6700    0.1900],'LineWidth',1);
plot(ax_map(4),[x_sq,fliplr(x_sq),x_sq(1)],[y_sq(1),y_sq(1),y_sq(2),y_sq(2),y_sq(1)],'Color',[0.4700    0.6700    0.1900],'LineWidth',1);

scatter(ax_map(2),stim_x(2),stim_x(1),80,[1 0 0],'p','filled');

im_map = imagesc(ax_map(4),[0,1],[0,1],O.param.f(reshape(mean(R.Var.V(:,3000:3100) - R.Var.phi(:,3000:3100),2),O.n)));
uistack(im_map,'bottom')
ax_map(4).XLim = ax_map(2).XLim;
ax_map(4).YLim = ax_map(2).YLim;
caxis(ax_map(4),[0 0.06]);
colormap(ax_map(4),[1 1 1;1 0 0]);

%%

scatter(ax_map(4),stim_x(2),stim_x(1),80,[1 0 0],'p','filled');

%% Import pre-made schematics
% ImPatient = imread('D:\Paper Drafts\Model paper\Figures\1BPatient','png','BackgroundColor',[1 1 1]);
% image(ax_map(1),ImPatient);
% set(ax_map(1),'XLim',[700 1500],'YLim',[700 1500]);
% ImPatient = imread('D:\Paper Drafts\Model paper\Figures\1CPatient','png','BackgroundColor',[1 1 1]);
% image(ax_map(3),ImPatient);
% set(ax_map(3),'XLim',[700 1500],'YLim',[700 1500]);
%% After you finish
Char = 'ABC';
for i = 1:numel(Char)
    An(i) = annotation('TextBox',[0.1 0.9 0.1 0.1],'Unit','normalized','String',Char(i),'FontSize',14,'FontWeight','bold','EdgeColor','none');
end
An(2).String = 'B    Slow Time Scale';
An(3).String = 'C    Fast Time Scale';

%% Make a video 
VIDEO = VideoWriter('Figure1Video');
VIDEO.Quality  = 55;
VIDEO.FrameRate = 100;
open(VIDEO)
f = figure('Color',[1 1 1]);f.Units = 'inches';f.Position = [1 1 1.9 1.5];
ax_video = axes('Position',[0 0.05 1 0.9]);
rate = (R.Var.V - R.Var.phi);
ImData = O.param.f(reshape(rate(:,1),100,100))/0.2;
ImData(ImData == 0) = inf;
Im = imagesc(ax_video,ImData);
ax_video.XAxis.Visible = 'off';
ax_video.YAxis.Visible = 'off';
ax_video.Units = 'inches';
ax_video.Position(3) = ax_video.Position(4);
ax_video.Units = 'normalized';
tx_video = text(ax_video,103,91,{'0.00','sec'},'HorizontalAlignment','center');
caxis([0 0.5])
H = hot(360);
colormap(ax_video,H);
cbar_video = colorbar('peer',ax_video,'east');
cbar_video.Position = [0.84 0.05 0.05 0.9];
cbar_video.Ticks = [0 0.5];
cbar_video.TickLabels{2} = '.5';
yo = ylabel(cbar_video,'f_m_a_x');
yo.Position(1) = 1;

for i = 1:8200
    ImData = O.param.f(reshape(rate(:,i),100,100))/0.2;
    ImData(ImData == 0) = inf;
    Im.CData = ImData;
    drawnow;
    im = getframe(gcf);
    tx_video.String = {num2str(i/100,'%.2f'),'sec'};
    writeVideo(VIDEO,im);
end
close(VIDEO);