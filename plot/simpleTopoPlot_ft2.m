function simpleTopoPlot_ft2(data, layoutFileName,labelsFlag,cmap,newF,statF)
%simpleTopoPlot_ft(data, layoutFileName, labelsFlag,cmap)
% in this version, data is a matrix cond*elec
if nargin<3
    labelsFlag='off';
end
if nargin<4
    cmap='jet';
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
dat=[];
dat.avg=mean(data,1);
dat.dimord='chan';
dat.label=layout.label(1:size(data,2));
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

if statF(1)
  
end