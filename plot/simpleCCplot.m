function simpleCCplot(MatCC,ChanLabels,lay,limPlot,newP)
%(MatCC,ChanLabels,lay,limPlot,newP)
if nargin<5 || newP==1
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

%%%%%%%%%%%%%%%%
if ndims(MatCC)==2
    MatCCnorm=floor((MatCC-min(min(MatCC)))/max(max(MatCC-min(min(MatCC))))*63)+1;
    cmap=colormap(jet);
    for i=1:length(ChanLabels)
        for j=i+1:length(ChanLabels)
            xbeg = Coor2D(1,i);
            ybeg = Coor2D(2,i);
            xend = Coor2D(1,j);
            yend = Coor2D(2,j);
            
            x = [xbeg xend]';
            y = [ybeg yend]';
            % h = line(x, y);
            if limPlot(1)==2
                if limPlot(3)==0
                    plotF=abs(MatCC(i,j))>limPlot(2);
                elseif limPlot(3)==1
                    plotF=(MatCC(i,j))>limPlot(2);
                elseif limPlot(3)==-1
                    plotF=(MatCC(i,j))<limPlot(2);
                end
            elseif limPlot(1)==3
                temp=(MatCC);
                temp=reshape(temp,1,numel(temp));
                temp(isnan(temp))=[];
                [order idx]=sort(temp,'descend');
                thrToKeep=order(limPlot(2));
                plotF=(MatCC(i,j))>=thrToKeep;
            else
                plotF=1;
            end
            if plotF
                line(x, y,'Color', cmap(MatCCnorm(i,j),:));
            end
            %                 h = patch(x, y, 1);
            %
            %         set(h, 'LineWidth', MatCCnorm(i,j));
        end
    end
    
elseif ndims(MatCC)==3
    meanMatCC=squeeze(mean(MatCC,1));
    [h sigMatCC]=ttest(MatCC,0,[],[],1);
    sigMatCC=squeeze(sigMatCC);
    MatCCnorm=floor((meanMatCC-min(min(meanMatCC)))/max(max(meanMatCC-min(min(meanMatCC))))*63)+1;
    cmap=colormap(jet);
    for i=1:length(ChanLabels)
        for j=i+1:length(ChanLabels)
            if limPlot(1)==1
                plotF=sigMatCC(i,j)<limPlot(2);
            elseif limPlot(1)==2
                if limPlot(3)==0
                plotF=abs(meanMatCC(i,j))>limPlot(2);
                elseif limPlot(3)==1
                plotF=(meanMatCC(i,j))>limPlot(2);
                elseif limPlot(3)==-1
                plotF=(meanMatCC(i,j))<limPlot(2);
                end
                            elseif limPlot(1)==3
                temp=(meanMatCC);
                temp=reshape(temp,1,numel(temp));
                                temp(isnan(temp))=[];
[order idx]=sort(temp,'descend');
                thrToKeep=order(limPlot(2));
                plotF=(meanMatCC(i,j))>=thrToKeep;
                else
                plotF=1;
            end
            if plotF
                xbeg = Coor2D(1,i);
                ybeg = Coor2D(2,i);
                xend = Coor2D(1,j);
                yend = Coor2D(2,j);
                
                x = [xbeg xend]';
                y = [ybeg yend]';
                % h = line(x, y);
                line(x, y,'Color', cmap(MatCCnorm(i,j),:));
                
                %                 h = patch(x, y, 1);
                %
                %         set(h, 'LineWidth', MatCCnorm(i,j));
            end
        end
    end
end

axis('equal')
xlim([-0.6 0.6])
ylim([-0.6 0.6])
set(gcf,'Color','w')
set(gca,'Xcolor','w','Ycolor','w')
box off