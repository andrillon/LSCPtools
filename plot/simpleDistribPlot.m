function [xout nout]=simpleDistribPlot(x,bins,Color,newF,smF)
% [xout nout]=simpleDistribPlot(x,bins,Color,newF,smF)

if nargin<5
    smF=0;
end
if nargin<4
    newF=1;
end
if nargin<3
    Color='b';
end
if nargin<2
    bins=10;
end

if newF
    figure;
end
hold on;

if isvector(x) && ~iscell(x)
    if ~isempty(bins)
        if length(bins)==1
            [nout, xout]=hist(x,bins);
            nout=nout/sum(nout)*100;
        else
            [nout, bin]=histc(x,bins);
            nout=nout/sum(nout)*100;
            xout=bins;
        end
    else
        [nout, xout]=hist(x,10);
        nout=nout/sum(nout)*100;
    end
    if smF==0
        plot(xout,nout,'Color',Color,'LineWidth',2);
    else
        plot(xout,fastsmooth(nout,smF,3,0),'Color',Color,'LineWidth',2);
    end
elseif ismatrix(x) && ~iscell(x)
    for nline=1:size(x,1)
        if ~isempty(bins)
            if length(bins)==1
                [nout(nline,:), xout(nline,:)]=hist(x(nline,:),bins);
            else
                [nout(nline,:), bin]=histc(x(nline,:),bins);
                xout(nline,:)=bins;
            end
        else
            [nout(nline,:), xout(nline,:)]=hist(x(nline,:),10);
        end
        nout(nline,:)=nout(nline,:)/sum(nout(nline,:))*100;
    end
    for nline=1:size(x,1)
        if smF==0
            plot(xout(nline,:),nout(nline,:),'Color',Color{nline},'LineWidth',2);
        else
            plot(xout(nline,:),fastsmooth(nout(nline,:),smF,3,0),'Color',Color{nline},'LineWidth',2);
        end
    end
elseif iscell(x)
    for nline=1:length(x)
        if ~isempty(bins)
            if length(bins)==1
                [nout{nline}, xout{nline}]=hist(x{nline},bins);
            else
                [nout{nline}, bin]=histc(x{nline},bins);
                xout{nline}=bins;
            end
        else
            [nout{nline}, xout{nline}]=hist(x(nline,:),10);
        end
        nout{nline}=nout{nline}/sum(nout{nline})*100;
        if smF==0
            plot(xout{nline},nout{nline},'Color',Color{nline},'LineWidth',2);
        else
            plot(xout{nline},fastsmooth(nout{nline},smF,3,0),'Color',Color{nline},'LineWidth',2);
        end
    end
end
