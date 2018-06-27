function [outPos]=simpleTopoPlot2(data, pos, labels,labelsFlag,cmap,newF,lay,baryF)
%simpleTopoPlot(data, pos, labels,labelsFlag,cmap,newF,lay)

if nargin<4
    labelsFlag=1;
end
if nargin<5  || isempty(cmap)
    cmap='jet';
end
if nargin<6 || newF==1
    figure;
end
if nargin<7 || isempty(lay)
    lay.outline=[];
end
if nargin<8
    baryF=0;
end
hold on;

clim    = [-nanmax(abs(data(:)))-( nanmax(data(:))-nanmin(data(:)) )/63 , nanmax(abs(data(:)))];

colormap(cmap)
if size(data,1)>size(data,2)
    data=data';
end

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
dataI      = griddata(cpos(1,:)',cpos(2,:)',full(double(cdata')),XI,YI,'cubic');


try
    imagescnan((dataI))
catch
    imagesc((dataI))
end
hold on;
contour((dataI),...
    'linecolor',0.5.*ones(3,1))

if ~isempty(lay.outline)
    for i=1:length(lay.outline)
        if ~isempty(lay.outline{i})
            X = lay.outline{i}(:,1)+0.5;
            Y = lay.outline{i}(:,2)+0.5;
            hold on
            line(100*X, 100*Y, 'color', 'k', 'linewidth', 3, 'linestyle', '-');
        end
    end
end
fpos = cpos;
fpos(1,:) = fpos(1,:) - xmin;
fpos(2,:) = fpos(2,:) - ymin;
fpos(1,:) = fpos(1,:)./(dx);
fpos(2,:) = fpos(2,:)./(dy);
% fpos(2,:) = 100-fpos(2,:); 

if labelsFlag>0
    plot(fpos(1,:),fpos(2,:),...
        'ko');
    if labelsFlag==2
    text(fpos(1,:),fpos(2,:),clabels);
    end
end
outPos=fpos;

if baryF
%     cminmax=(cdata-min(cdata))./max(cdata-min(cdata))-0.5;
%     bary=nanmean(100*cpos.*repmat(cminmax,2,1),2);
cY = mean(Y(M==1));
cX = mean(X(M==1));

    scatter(bary(1),bary(2),'Marker','+','MarkerEdgeColor','k','SizeData',72)
end
% xlim([-10 110])
% ylim([-10 110])
axis('equal')
caxis([-1 1]*nanmax(abs(data)))

set(gcf,'Color','w')
set(gca,'Xcolor','w','Ycolor','w')