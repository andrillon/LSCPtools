function simpleCCplot2(MatCC,ChanLabels,lay,newP)
%(MatCC,ChanLabels,lay,limPlot,newP)
if nargin<4 || newP==1
    figure;
end

% plot head template
Coor2D=lay.pos(match_str(ChanLabels,lay.label),:)';
scatter(Coor2D(1,:), Coor2D(2,:),'k.');

for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
        X = lay.outline{i}(:,1);
        Y = lay.outline{i}(:,2);
        line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
    end
end

%
MatCCnorm=floor((MatCC-min(min(MatCC)))/max(max(MatCC-min(min(MatCC))))*63)+1;
cmap=colormap('jet');
for i=1:length(ChanLabels)
    for j=i+1:length(ChanLabels)
        xbeg = Coor2D(1,i);
        ybeg = Coor2D(2,i);
        xend = Coor2D(1,j);
        yend = Coor2D(2,j);
        
        x = [xbeg xend]';
        y = [ybeg yend]';
        line(x, y,'Color', cmap(MatCCnorm(i,j),:));
        
    end
end
axis('equal')
xlim([-0.6 0.6])
ylim([-0.6 0.6])
set(gcf,'Color','w')
set(gca,'Xcolor','w','Ycolor','w')
h=colorbar;
set(h,'Ytick',[0 1],'YTicklabel',[min(min(MatCC)) max(max(MatCC))]);
box off