function [databin]=simpleHistPlot(data,bins,type,Prop,normF,newF)
% simpleHistPlot(data,bins,type,Prop)
% data 1*n or m*n
% bins (number or bin range)
% line or bar
%   line: colo line, linewidth, linestyle, smoothflag
%   bar: [colorface; coloredge], linewidth, barwidth
% normF: normalize or not
% newF: naw figure

if exist('normF')==0
    normF=0;
end
if exist('newF')==0
    newF=1;
end

if newF
    figure;
else
    hold on;
end
set(gcf,'Color','w')
set(gca,'FontSize',16,'FontWeight','bold')

if isempty(type)
    type='bar';
end
if isempty(type)
    if strcmp(type,'bar')
        Prop={['k';'k'],2,1};
    elseif strcmp(type,'line')
        Prop={['k';'k'],2,'-',0};
    end
end
if size(Prop{1},1)==1
    Prop{1}(2,:)=Prop{1}(1,:);
end

if min(size(data))==1
    %%% hist
    if numel(bins)
        [n, outbin]=hist(data,bins);
    else
        [n, outbin]=histc(data,bins);
    end
    if normF
        n=n/sum(n)*100;
    end
    
    %%% plot
    if strcmp(type,'bar')
        bar(outbin,n,'FaceColor',Prop{1}(1,:),'EdgeColor',Prop{1}(2,:),'LineWidth',Prop{2},'BarWidth',Prop{3})
        databin=n;
    elseif strcmp(type,'line')
        databin(1,:)=n;
        if Prop{4}~=0
            n=fastsmooth(n,Prop{4},3,0);
        end
        databin(2,:)=n;
        plot(outbin,n,'Color',Prop{1}(1,:),'LineWidth',Prop{2},'LineStyle',Prop{3})
    end
else
    for nL=1:size(data,1)
        %%% hist
        if numel(bins)
            [n(nL,:), outbin]=hist(data(nL,:),bins);
        else
            [n(nL,:), outbin]=histc(data(nL,:),bins);
        end
        if normF
            n(nL,:)=n(nL,:)/sum(n(nL,:))*100;
        end
    end
    %%% plot
    if strcmp(type,'bar')
        meann=nanmean(n);
        bar(outbin,meann,'FaceColorBar',Prop{1}(1,:),'EdgeColorBar',Prop{1}(2,:),'LineWidth',Prop{2},'BarWidth',Prop{3})
        errorbar(outbin,meann,nansem(n),'Color',Prop{1}(1,:),'LineWidth',Prop{2})
        databin=meann;
    elseif strcmp(type,'line')
        meann=nanmean(n);
        databin(1,:)=meann;
        if Prop{4}~=0
            meann=fastsmooth(meann,Prop{4},3,0);
            databin(2,:)=meann;
        end
        plot(outbin,meann,'Color',Prop{1}(1,:),'LineWidth',Prop{2},'LineStyle',Prop{3})
        errorbar(outbin,nanmean(n),nansem(n),'Color',Prop{1}(1,:),'LineWidth',Prop{2})
    end
end

if normF
    ylabel('%')
else
    ylabel('#')
end