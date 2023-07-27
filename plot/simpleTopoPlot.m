simplfunction simpleTopoPlot(data, pos, labels,labelsFlag,cmap)
%simpleTopoPlot(data, pos, labels,labelsFlag,cmap)

if nargin<4
    labelsFlag=1;
end
if nargin<5
    cmap=[];
end
ParentAxes = [];
f = [];
% clim    = [-max(abs(data(:)))-( max(data(:))-min(data(:)) )/63 , max(abs(data(:)))];
figName = 'Image Scalp data';
noButtons = 0;


% exclude channels ?
goodChannels = find(~isnan(pos(1,:)));
pos          = pos(:,goodChannels);
data            = data(goodChannels,:);
labels    = labels(goodChannels);

cdata         = data;
cpos       = pos;
clabels = labels;

xmin    = min(cpos(1,:));
xmax    = max(cpos(1,:));
dx      = (xmax-xmin)./100;
ymin    = min(cpos(2,:));
ymax    = max(cpos(2,:));
dy      = (ymax-ymin)./100;
x       = xmin:dx:xmax;
y       = ymin:dy:ymax;
[XI,YI] = meshgrid(x,y);
dataI      = griddata(cpos(1,:)',cpos(2,:)',full(double(cdata')),XI,YI);


f=figure(...
    'name',figName,...
    'color',[1 1 1]);%,...
%     'deleteFcn',@dFcn);
ParentAxes = axes('parent',f);

COLOR = get(f,'color');
d.hi = imagescnan(flipud(dataI),...
    'CDataMapping','scaled',...
    'Parent',ParentAxes);
set(ParentAxes,'nextPlot','add',...
    'tag','spm_eeg_plotScalpData')
if length(unique(dataI)) ~= 1
    [C,d.hc] = contour(ParentAxes,flipud(dataI),...
        'linecolor',0.5.*ones(3,1));
end
caxis(ParentAxes,clim);
if isempty(cmap)
col = jet;
else
    col=cmap;
end
% col(1,:) = COLOR;
colormap(ParentAxes,col)

d.cbar = colorbar('peer',ParentAxes);

axis(ParentAxes,'off')
axis(ParentAxes,'equal')
axis(ParentAxes,'tight')

fpos = cpos;
fpos(1,:) = fpos(1,:) - xmin;
fpos(2,:) = fpos(2,:) - ymin;
fpos(1,:) = fpos(1,:)./(dx);
fpos(2,:) = fpos(2,:)./(dy);
fpos(2,:) = 100-fpos(2,:);  % for display purposes (flipud imagesc)

fpos = cpos;
fpos(1,:) = fpos(1,:) - xmin;
fpos(2,:) = fpos(2,:) - ymin;
fpos(1,:) = fpos(1,:)./(dx);
fpos(2,:) = fpos(2,:)./(dy);
fpos(2,:) = 100-fpos(2,:);  % for display purposes (flipud imagesc)

figure(f);
if labelsFlag
d.hp = plot(ParentAxes,...
    fpos(1,:),fpos(2,:),...
    'ko');
    d.ht = text(fpos(1,:),fpos(2,:),clabels,...
        'Parent',ParentAxes,...
        'visible','on');
end
axis(ParentAxes,'image')

