function [pV,H,H2]=simpleSpreadPlot(Pos,data,colorBar,widthBar,sigFlag,widthLine)
% hbar=simpleBarPlot(Pos,data,colorBar,widthBar,colorError,sigFlag)
pV=[];
if nargin<6 || isempty(sigFlag)
    sigFlag={0 0 0};
end
if nargin<7
    widthLine=3;
end
if size(colorBar,1)==1
    if colorBar(1)=='w'
        colorBar(2,:)='k';
    else
        colorBar(2,:)=colorBar(1,:);
    end
end
hold on;
if sigFlag{1}==1
    
    H = boxplot(data,'color',colorBar(1,:),'positions',Pos,'widths',widthBar, 'symbol','');%,'color',Colors2([1 3 2 4]))
    H2= plotSpread(data,[],[],[],[],Pos);
    set(H2{1},'LineWidth',2,'Marker', 'o','Color',colorBar(1,:),'MarkerFaceColor',colorBar(2,:))
    set(H,'LineWidth',widthLine)
    
    if length(data)==length(sigFlag{2}) || length(sigFlag{2})==1
        [h pV]=ttest(data,sigFlag{2});
        fprintf('... paired t-test p=%1.5f\n',pV);
    else
        [h pV]=ttest2(data,sigFlag{2});
        fprintf('... non-paired t-test p=%1.5f\n',pV);
    end
    if pV<sigFlag{3}
        if nanmean(data)>0
            plot(Pos,1.1*nanmax(data),'*k')
        else
            plot(Pos,-1.1*nanmax(data),'*k')
        end
    end
elseif sigFlag{1}==2
    
    H = boxplot(data,'color',colorBar(1,:),'positions',Pos,'widths',widthBar, 'symbol','');%,'color',Colors2([1 3 2 4]))
    H2= plotSpread(data,[],[],[],[],Pos);
    set(H2{1},'LineWidth',2,'Marker', 'o','Color',colorBar(1,:),'MarkerFaceColor',colorBar(2,:))
    set(H,'LineWidth',widthLine)
    
    if length(data)==length(sigFlag{2}) || length(sigFlag{2})==1
        [pV h]=signrank(data,sigFlag{2});
        fprintf('... paired u-test p=%1.5f\n',pV);
    else
        [pV h]=ranksum(data,sigFlag{2});
        fprintf('... non-paired u-test p=%1.5f\n',pV);
    end
    if pV<sigFlag{3}
         if nanmean(data)>0
            plot(Pos,1.1*nanmax(data),'*k')
        else
            plot(Pos,-1.1*nanmax(data),'*k')
        end
    end
else
    H = boxplot(data,'color',colorBar(1,:),'positions',Pos,'widths',widthBar, 'symbol','');%,'color',Colors2([1 3 2 4]))
    H2= plotSpread(data,[],[],[],[],Pos);
    set(H2{1},'LineWidth',2,'Marker', 'o','Color',colorBar(1,:),'MarkerFaceColor',colorBar(2,:))
    set(H,'LineWidth',widthLine)
    
end