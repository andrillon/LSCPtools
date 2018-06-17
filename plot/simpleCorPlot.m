function [stats hdl]=simpleCorPlot(X,Y,Prop,corType,newF)
% [rho pV hdl]=simpleCorPlot(X,Y,Prop,corType,newF)

rho=[]; pV=[];
if nargin<3 || isempty(Prop)
    Prop={'o','b','b',72};
end
if nargin<4
    corType=[];
end
if nargin<5
    newF=0;
end
% plot scatter
if newF
    figure;
end
format_fig;
% oldaxis=axis;
if length(Prop)==4
    if size(Prop{3},1)==1
        if isempty(Prop{2})
            hdl=scatter(X,Y,'LineWidth',3,'Marker',Prop{1},'MarkerEdgeColor',Prop{3},'SizeData',Prop{4},'LineWidth',4);
        else
            hdl=scatter(X,Y,'LineWidth',3,'Marker',Prop{1},'MarkerFaceColor',Prop{2},'MarkerEdgeColor',Prop{3},'SizeData',Prop{4});
        end
    else
        if isempty(Prop{2})
            hdl=scatter(X,Y,[],Prop{3},'LineWidth',3,'Marker',Prop{1},'SizeData',Prop{4},'LineWidth',4);
        else
            hdl=scatter(X,Y,[],Prop{3},'LineWidth',3,'Marker',Prop{1},'MarkerFaceColor',Prop{2},'SizeData',Prop{4});
        end
    end
elseif length(Prop)==5
    [Xorder,idxorder]=sort(X);
    Yorder=Y(idxorder);
    step=floor(length(Xorder)/Prop{5});
    for k=step:step:length(Xorder)
        meanx=nanmean(Xorder(k-step+1:k));
        semx=sem(Xorder(k-step+1:k));
        meany=nanmean(Yorder(k-step+1:k));
        semy=sem(Yorder(k-step+1:k));
        
        line([1 1]*meanx,[-1 1]*semy+meany,'Color',Prop{2},'LineWidth',2)
        line([-1 1]*semx+meanx,[1 1]*meany,'Color',Prop{2},'LineWidth',2)
    end
end
% axis(oldaxis);
hold on;
% correlation
if size(X,1)<size(X,2)
    X=X';
end
if size(Y,1)<size(Y,2)
    Y=Y';
end
if ~isempty(corType)
    if strcmp(corType,'robustfit')
        try
            [b,rstats] = robustfit(X,Y);
            plot(xlim,b(1)+b(2)*xlim,'Color',Prop{2},'LineStyle','--','LineWidth',3);
            fprintf('... ... robust fit coefficient: Y = %g * X + %g (%1.3f, %1.3f)\n',b(2),b(1),rstats.p(2),rstats.p(1))
            stats=[b rstats.p];
        catch
            stats=nan(2,2);
            fprintf('... ... could not fit\n')
        end
    elseif strcmp(corType,'pearson_nofit')
          [rho pV]=corr(X,Y,'type','pearson','rows','pairwise');
        stats=[rho pV];
        
        if newF
            title(sprintf('%s correlation: rho=%g p=%1.3f',corType,rho,pV))
        else
            title(sprintf('%s: r=%1.2f p=%1.3f',corType,rho,pV))
            fprintf('%s correlation: r=%1.2f p=%1.3f\n',corType,rho,pV)
        end
    else
        [rho pV]=corr(X,Y,'type',corType,'rows','pairwise');
        stats=[rho pV];
        
        if pV<0.1
           [b,rstats] = robustfit(X,Y);
            plot(xlim,b(1)+b(2)*xlim,'Color',Prop{3},'LineStyle','--','LineWidth',3);
        else
             [b,rstats] = robustfit(X,Y);
            plot(xlim,b(1)+b(2)*xlim,'Color',Prop{3},'LineStyle',':','LineWidth',3);
        end
        %    lin=linspace(min(X),max(X));
        %    hold on
        %    plot(lin,lin*rho,'k');
        if newF
            title(sprintf('%s correlation: r=%1.2f p=%1.3f',corType,rho,pV))
        else
            title(sprintf('%s: r=%1.2f p=%1.3f',corType,rho,pV))
            fprintf('%s correlation: r=%1.2f p=%1.3f\n',corType,rho,pV)
        end
    end
end
