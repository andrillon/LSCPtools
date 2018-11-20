function handles = simple_PiratePlot(X,Y,width,smoothingFactor,cutEdges,pcolor,spread,limits,scatterargs)
% Trace one or several custom Pirate Plot. Circles represent individual
% data, middle bar the mean of the distribution, shaded area the mean +/-
% the standard error of the mean (corrected for small samples; This is not
% the standard Pirate Plot, shaded area should represent the Bayesian
% Highest Density Interval (HDI)) and the lines the smoothed distribution
% of the data.
%
%     gyt_PiratePlot(Y)
%     gyt_PiratePlot(X,Y)
% Specify a unique vector 'Y' (one Pirate Plot) or a unique matrix 'Y' (as
% many Pirate Plot as the length of 'Y'). Alternatively, specify 'X' as a
% scalar (for one plot) or as a vector and 'Y' must then be a vector (for
% one plot) or a matrix with at least one dimension of equal size than 'X'.
%
%     gyt_PiratePlot(X,Y,width)
% 'width' argument specify the width of the pirate plot on the X axis of
% the plot (default: width = 1).
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor)
% 'smoothingFactor' is a scalar that control the smoothing of the
% distribution reprensented by the lines (default: smoothingFactor = 2).
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor,cutEdges)
% 'cutEdges' is a single character that control whether 0 values of the
% smoothed distribution should be displayed or not. 'cutEdges' can be
% either 'y' for yes (cut edges), 'p' for partial (cut only at the ends of
% the distribution) or 'n' for no (no cut) (default: cutEdges = 'y').
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor,cutEdges,color)
% 'color' is a vector of 3 elements or a matrix of Nx3 elements with N the
% number of elements in X (default: color = [0,0,1]).
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor,cutEdges,color,spread)
% 'spread' is a single character that control whether individual data
% should be shown. 'spread' can either be 'y' for yes or 'n' for no
% (default: spread = 'y').
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor,cutEdges,color,spread,limits)
% 'limits' is a 2-elements vector that control the boundaries of the
% smoothed distribution. If 'limits' is empty, there are no boundaries to
% the smoothed distribution (default: limits = []).
%
%     gyt_PiratePlot(X,Y,width,smoothingFactor,cutEdges,color,spread,limits,scatterargs)
% 'scatterargs' is a cell array controling the graphical properties of the
% scatter plot used to show individual data (default: scatterargs = {}).
%
%
% 
% This script was created from the guidelines of Nathaniel D. Phillips and
% his R package "yarrr", all detailed in the excellent article:
% https://www.r-bloggers.com/the-pirate-plot-an-r-pirates-favorite-plot/
% 
% Author : Guillaume Legendre  -- 31/07/2017 --
% UNIGE




if nargin<2
    Y = X;
    if sum(size(Y)>1)==1
        X = 1;
    else
        X = 1:size(Y,1);
    end
end
if nargin<3 || isempty(width)
    width = 1;
end
if nargin<4 || isempty(smoothingFactor)
    smoothingFactor = 2;
end
if nargin<5 || isempty(cutEdges)
    cutEdges = 'y';
end
if nargin<6 || isempty(pcolor)
    pcolor = [0,0,1];
end
if nargin<7 || isempty(spread)
    spread = 'y';
end
if nargin<8 || isempty(limits)
    limits = [];
end
if nargin<9 || isempty(scatterargs)
    scatterargs = [];
end

if ndims(Y)>2
    error('Cannot work with matrices of more than 2 dimensions')
end
if ~any(numel(X)==size(Y))
    error('X must contain as many elements than one of the dimension of Y')
end
if size(pcolor,1)~=1 && size(pcolor,1)~=numel(X)
    error('When specifying colors, ''color'' should be a vector with 3 elements or a matrix of size Nx3 with N the number of elements of X.')
end


if any(find(size(Y)==numel(X))==1)
    Y = Y';
end


hold on


for iP = 1:numel(X)
    
    yplot = Y(:,iP)';
    xplot = X(iP);
    if size(pcolor,1)==numel(X)
        ccolor = pcolor(iP,:);
    else
        ccolor = pcolor;
    end
    
    darkcolor = ccolor-[0.1,0.1,0.1];
    darkcolor(darkcolor<0) = 0;
    brightcolor = ccolor+[0.35,0.35,0.35];
    brightcolor(brightcolor>1) = 1;
    
    divisions = 10;
    
    binsize = iqr(yplot)/(sum(~isnan(yplot)).^(1/3))/2;
    if binsize==0
        binsize=0.02;
    end
    totbin = ceil((nanmax(yplot)-nanmin(yplot))/binsize);
    totbin = totbin*2;
    
    allbins = (nanmedian(yplot)-totbin*binsize/2)+(0:totbin)*binsize;
    freqs = histc(yplot,allbins);
    freqs = freqs(1:(end-1))/numel(yplot);
    x = repmat(freqs,divisions,1);
    x = x(:)';
    
    N = round(divisions/smoothingFactor);
    defaultKernel = ones(1,N)/N;
    kernel = [zeros(1,N),defaultKernel,zeros(1,N)];
    kernel = conv(kernel,defaultKernel,'same');
    kernel = conv(kernel,defaultKernel,'same');
    xx = conv(x,kernel,'same');
    xx = xx/max(xx)*width;



    allpoints = linspace(nanmedian(yplot)-totbin(end)*binsize/2,nanmedian(yplot)+totbin(end)*binsize/2,totbin*divisions);
    if ~isempty(limits)
        inds = find((allpoints>=limits(1)).*(allpoints<=limits(2)));
        xx = xx(inds);
        allpoints = allpoints(inds);
    end
    if strcmpi(cutEdges,'y')
        divs = cumsum([1,diff(xx>0)~=0]).*(xx>0);
        divsI = setdiff(unique(divs),0);
        handles = nan(1,numel(divsI)+3);
        for nC = 1:numel(divsI)
            handles(nC) = fill(xplot+[-xx(divs==divsI(nC)),fliplr(xx(divs==divsI(nC)))],[allpoints(divs==divsI(nC)),fliplr(allpoints(divs==divsI(nC)))],'w','edgecolor',darkcolor);
        end
    elseif strcmpi(cutEdges,'n')
        handles = nan(1,4);
        handles(1) = fill(xplot+[-xx,fliplr(xx)],[allpoints,fliplr(allpoints)],'w','edgecolor',darkcolor);
    elseif strcmpi(cutEdges,'p')
        divs = cumsum([1,diff(xx>0)~=0]).*(xx<=0);
        divsI = setdiff(unique(divs),0);
        keptPart = ~ismember(divs,divsI([1,end]));
        handles = nan(1,4);
        handles(1) = fill(xplot+[-xx(keptPart),fliplr(xx(keptPart))],[allpoints(keptPart),fliplr(allpoints(keptPart))],'w','edgecolor',darkcolor);
    end

    totInt = trapz(xx);
    [sortxx] = sort(xx,'descend');
    curralpha = 0;
    inc = 2;
    while curralpha<0.95
        tempxx = xx;
        tempxx((xx-sortxx(inc))<=0)=0;
        tempInt = trapz(tempxx);
        curralpha = tempInt/totInt;
        inc = inc+1;
    end

    [~,Ind] = nanmin(abs(allpoints-nanmean(yplot)));
%     SEMsize = tinv(0.975,sum(~isnan(yplot))-1)*(nanstd(yplot)/sqrt(sum(~isnan(yplot))-1));
    SEMsize = nanstd(yplot)/sqrt(sum(~isnan(yplot))-1);
    firstI = find(diff((nanmean(yplot)-SEMsize)>allpoints)~=0);
    if isempty(firstI)
        firstI=1;
    end
    lastI = find(diff((nanmean(yplot)+SEMsize)>allpoints)~=0);
    if isempty(lastI)
        lastI=numel(allpoints)-1;
    end
    firstP = ((nanmean(yplot)-SEMsize)-allpoints(firstI))/(allpoints(firstI+1)-allpoints(firstI))*(xx(firstI+1)-xx(firstI))+xx(firstI);
    lastP = ((nanmean(yplot)+SEMsize)-allpoints(lastI))/(allpoints(lastI+1)-allpoints(lastI))*(xx(lastI+1)-xx(lastI))+xx(lastI);
    SEMxplotspan = [firstP,xx((firstI+1):lastI),lastP];
    SEMyplotspan = [nanmean(yplot)-SEMsize,allpoints((firstI+1):lastI),nanmean(yplot)+SEMsize];
    
       handles(end-2) = fill(xplot+[-SEMxplotspan,fliplr(SEMxplotspan)],[SEMyplotspan,fliplr(SEMyplotspan)],brightcolor,'edgecolor','none');
 if strcmpi(spread,'y')
        
        [yplot,Iyy] = sort(yplot);

        rap = max(diff(yplot));

        yy = [];
        minidist = rap;

        for iX=1:numel(yplot)
            if ~isempty(yy)
                closepoints = abs(yy(1,:)-yplot(iX))<minidist;
            end
            if isempty(yy) || ~any(closepoints)
                pointCoord = [yplot(iX);0];
            else
                yadd = sqrt((((minidist.^2)*ones(1,size(yy,2)))-(yy(1,:)-yplot(iX)).^2).*closepoints);
                minYs = yy(2,:)-yadd;
                maxYs = yy(2,:)+yadd;
                minYs = minYs(closepoints==1);
                maxYs = maxYs(closepoints==1);
                divX = [minYs,maxYs];
                matX = repmat(divX,[sum(closepoints),1]);
                divXMark = divX(all(((matX-repmat(minYs',[1,numel(divX)]))<=0)+((matX-repmat(maxYs',[1,numel(divX)]))>=0)));
                [~,iY] = min(abs(divXMark-0));
                tempy = divXMark(iY);
                tempy = tempy+normrnd(0,0.25); %*minidist;
                pointCoord = [yplot(iX);tempy];
            end
            yy = [yy,pointCoord];
        end
        [~,Iyy2] = sort(Iyy);
        yy = yy(2,Iyy2);
        yy = yy/max(yy)*width*0.5;
        
        handles(end) = scatter(xplot+yy,yplot,2,'markerfacecolor',ccolor,'markeredgecolor','k');
        if ~isempty(scatterargs)
            set(handles(end),scatterargs{:});
        end
    end
    handles(end-1) = plot(xplot+[-xx(Ind),xx(Ind)],[nanmean(yplot),nanmean(yplot)],'color',darkcolor,'linewidth',1);

end


