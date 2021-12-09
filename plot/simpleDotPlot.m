function [hdot pV]=simpleDotPlot(Pos,data,sizeDot,colorBar,widthBar,colorError,markerType,sigFlag,widthLine,spreadFlag,boxFlag,splitFlag)
% INPUT:
% - Pos         : x-position of the center of the bar
% - data        : vector of data to plot
% - colorBar    : color of the bar (string or RGB value or 2 RGB values
% (line and center)
% - widthBar    : width of the bar (e.g. 1)
% - colorError  : color of the error bar (e.g. 'r' for red)
% - sigFlag     : 3-elements cell
%                       0, 1 or 2: no stat, param stat or non-param stat
%                       H0 (one value or vector to compare)
%                       p-value threshold (e.g. 0.05)
% - widthLine   : witdht of the bar line (eg 1)

pV=[];
if nargin<8 || isempty(sigFlag)
    sigFlag={0};
end
if nargin<9
    widthLine=3;
end
if nargin<10
    spreadFlag=0;
end
if nargin<11
    boxFlag=0;
end
if nargin<12
    splitFlag=0;
end
if nargin<7
    markerType='o';
end
if size(colorBar,1)==1
    if colorBar(1)=='w'
        colorBar(2,:)='k';
    else
        colorBar(2,:)=colorBar(1,:);
    end
end
if isempty(sizeDot)
    sizeDot=36;
end
hold on;
if boxFlag
    line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1].*prctile(data,25),'Color',colorBar(1,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1].*prctile(data,75),'Color',colorBar(1,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos-0.1*widthBar]+splitFlag,[prctile(data,25) prctile(data,75)],'Color',colorBar(1,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1].*nanmean(data),'Color',colorBar(1,:),'LineWidth',widthLine+1)
    line([Pos+0.1*widthBar Pos+0.1*widthBar]+splitFlag,[prctile(data,25) prctile(data,75)],'Color',colorBar(1,:),'LineWidth',widthLine)
    
    patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]'+splitFlag,...
        [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
    %    minBoxPlot=prctile(data,25)-1.5*(prctile(data,75)-prctile(data,25));
    %    maxBoxPlot=prctile(data,75)+1.5*(prctile(data,75)-prctile(data,25));
    %    line([1 1]*Pos,[minBoxPlot prctile(data,25)],'Color',colorBar(1,:),'LineWidth',widthLine-1)
    %    line([1 1]*Pos,[prctile(data,75) maxBoxPlot],'Color',colorBar(1,:),'LineWidth',widthLine-1)
    %
    %    line([Pos-0.1*widthBar/2 Pos+0.1*widthBar/2],[1 1].*minBoxPlot,'Color',colorBar(1,:),'LineWidth',widthLine)
    %    line([Pos-0.1*widthBar/2 Pos+0.1*widthBar/2],[1 1].*maxBoxPlot,'Color',colorBar(1,:),'LineWidth',widthLine)
    if spreadFlag
        xspread=(rand(1,length(data))-0.5)*widthBar/4+Pos-splitFlag;
        hdot.ind=scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
    
else
    if spreadFlag
        xspread=(rand(1,length(data))-0.5)*widthBar/4+Pos-splitFlag;
        hdot.ind=scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
    
    if sigFlag{1}==1
        
        hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
        hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
        hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
        if length(data)==length(sigFlag{2}) || length(sigFlag{2})==1
            [h, pV, ~, stats]=ttest(data,sigFlag{2});
            fprintf('... paired t-test p=%1.5f (t(%g)=%2.3f)\n',pV,stats.df,stats.tstat);
        else
            [h, pV, ~, stats]=ttest2(data,sigFlag{2});
            fprintf('... non-paired t-test p=%1.5f (t(%g)=%2.3f)\n',pV,stats.df,stats.tstat);
        end
        if pV<sigFlag{3}
            if nanmean(data)>0
                plot(Pos+splitFlag,1.3*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'*k')
            else
                plot(Pos+splitFlag,-1.3*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'*k')
            end
        end
        hdot.mean=scatter(Pos,nanmean(data),'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorError,'LineWidth',widthLine,'SizeData',sizeDot);
    elseif sigFlag{1}==2
        hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError);
        hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError);
        hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError);
        
        if length(data)==length(sigFlag{2}) || length(sigFlag{2})==1
            [pV h]=signrank(data,sigFlag{2});
            fprintf('... paired u-test p=%1.5f\n',pV);
        else
            [pV h]=ranksum(data,sigFlag{2});
            fprintf('... non-paired u-test p=%1.5f\n',pV);
        end
        if pV<sigFlag{3}
            if nanmedian(data)>0
                plot(Pos+splitFlag,1.3*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'*k')
            else
                plot(Pos+splitFlag,-1.3*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'*k')
            end
        end
        hdot.mean=scatter(Pos+splitFlag,nanmedian(data),'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorError,'LineWidth',widthLine,'SizeData',sizeDot);
    elseif sigFlag{1}==0 && length(sigFlag)==3
        if length(data)==length(sigFlag{2}) || length(sigFlag{2})==1
            [h, pV, ~, stats]=ttest(data,sigFlag{2});
            fprintf('... paired t-test p=%1.5f (t(%g)=%2.3f)\n',pV,stats.df,stats.tstat);
        else
            [h, pV, ~, stats]=ttest2(data,sigFlag{2});
            fprintf('... non-paired t-test p=%1.5f (t(%g)=%2.3f)\n',pV,stats.df,stats.tstat);
        end
        if widthLine>1
            hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            
        else
            hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            
        end
        hdot.mean=scatter(Pos+splitFlag,nanmean(data),'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorError,'LineWidth',widthLine,'SizeData',sizeDot);
    elseif sigFlag{1}==4
        line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmedian(data),'LineWidth',widthLine-1,'Color',colorError)
        
        hdot=scatter(Pos+splitFlag,nanmean(data),'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorError,'LineWidth',widthLine,'SizeData',sizeDot);
    else
        if widthLine>1
            hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine-1,'Color',colorError);
            
        else
            hdot.line{1}=line([1 1]*Pos+splitFlag,[-1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            hdot.line{2}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[-1 -1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            hdot.line{3}=line([Pos-0.1*widthBar Pos+0.1*widthBar]+splitFlag,[1 1]*nanstd(data)/sqrt(sum(~isnan(data))-1)+nanmean(data),'LineWidth',widthLine,'Color',colorError);
            
        end
        hdot.mean=scatter(Pos+splitFlag,nanmean(data),'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorError,'LineWidth',widthLine,'SizeData',sizeDot);
    end
end