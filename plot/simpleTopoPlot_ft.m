function simpleTopoPlot_ft(data, layoutFileName,labelsFlag,cmap,newF,contourF)
%simpleTopoPlot_ft(data, layoutFileName, labelsFlag,cmap,newF)

if nargin<3
    labelsFlag='off';
end
if nargin<4 || isempty(cmap)
    cmap='parula';
end
if nargin<6
    contourF=1;
end
if nargin<5
    newF=1;
end
if newF
    figure;
end

if isstr(layoutFileName)
load(layoutFileName)
else
    layout=layoutFileName;
end
cfg = [];
cfg.layout = layout;
cfg.half=[];
% cfg.zlim = [-2.5 2.5]; 
cfg.comment='no';
cfg.marker=labelsFlag;
cfg.gridscale=256;
cfg.interpolation='v4';
if contourF==1
    cfg.style='both';
elseif contourF==1
    cfg.style='straight';
elseif contourF==2
    cfg.style='contour';
end
if newF
    cfg.figure='yes'
else
    cfg.figure='gcf';
end
dat=[];
dat.avg=data;
dat.dimord='chan';
dat.label=layout.label(1:length(data));
% if length(data)==65
%     cfg.channel=layout.label;
%     cfg.channel(match_str( cfg.channel,{'E62','E63'}))=[];
% dat.avg=data(setdiff(1:65,[62 63]));
% dat.label=layout.label(setdiff(1:65,[62 63]));
% 
% end
ft_topoplotER(cfg,dat);
colormap(cmap);

axis('equal')

set(gcf,'Color','w')
set(gca,'Xcolor','w','Ycolor','w')