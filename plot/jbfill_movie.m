function h=jbfill_movie(xdata,data,Color,Transp,speed,newF)
%%% jbfill_movie(data,Color,Transp,speed,newF)
%%%
%%% plot jbfill as a movie
%%% inputs:
%%% - data (repeats * time)
%%% - Color (color of your plot)
%%% - Transp (transpareny: 0 to 1)
%%% - speed (pause between plots in s)
%%% _ plot on a new figure or not
if newF==0
    figure(gcf);
    myaxis=axis;
else
    figure;
    xlim([xdata(1) xdata(end)])
    ylim([1.1*min(min(data)) 1.1*max(max(data))])
        myaxis=axis;

end
for n=1:length(xdata)
    cla;
    xlim(myaxis(1:2))
    ylim(myaxis(3:4))
    hold on;
    h=jbfill(xdata(1:n),mean(data(:,1:n))+std(data(:,1:n)),mean(data(:,1:n))-std(data(:,1:n)),Color,Color,0,Transp);
    plot(xdata(1:n),mean(data(:,1:n)),'LineWidth',2,'Color',Color)
    pause(speed)
end