function h=simpleTFplot(dat,fr,tp,corrplot,newplot)

% simpleTFplot(dat,fr,tp,corrplot,newplot)
if nargin<5
    newplot=0;
end
if newplot
    figure;
else
   figure(gcf);
end
format_fig;
if length(unique(diff(fr)))>2
    frtick=round(fr);
    frtickflag=1;
    fr=1:length(fr);
else
    frtickflag=0;
end
if isnumeric(dat)
    if corrplot
        datplot=log(dat./repmat(nanmean(dat(:,tp<0),2),1,size(dat,2)));
    else
        datplot=dat;
    end
    h=imagesc(tp,fr,datplot); set(gca,'YDir','normal');
    xlabel('Time')
    ylabel('Freq (Hz)')
    if newplot
        colorbar;
    end
    % clabel('Power (log)')
    caxis([-max(max(abs(datplot))) max(max(abs(datplot)))])
elseif iscell(dat)
    for n=1:length(dat)
        subplot(1,length(dat)+1,n)
        if corrplot
            datplot{n}=log(dat{n}./repmat(nanmean(dat{n}(:,tp<0),2),1,size(dat{n},2)));
        else
            datplot{n}=dat{n};
        end
        h=imagesc(tp,fr,datplot{n}); set(gca,'YDir','normal');
        xlabel('Time')
        ylabel('Freq (Hz)')
        if newplot
            colorbar;
        end
        % clabel('Power (log)')
        caxis([-max(max(abs(datplot{n}))) nanmax(nanmax(abs(datplot{n})))])
    end
    subplot(1,length(dat)+1,n+1)
    diffplot=datplot{2}-datplot{1};
    h=imagesc(tp,fr,diffplot); set(gca,'YDir','normal');
    xlabel('Time')
    ylabel('Freq (Hz)')
    if newplot
        colorbar;
    end
    % clabel('Power (log)')
    caxis([-max(max(abs(diffplot))) nanmax(nanmax(abs(diffplot)))])
end
if frtickflag
    downs=ceil(length(fr)/10);
    set(gca,'ytick',fr(1:downs:end),'yticklabel',frtick(1:downs:end))

end