function [rhopop, pVpop]=simpleCorPlotbin(X,Y,Prop,corType,newF)
% [rho pV hdl]=simpleCorPlot(X,Y,Prop,corType,newF)

rho=[]; pV=[];
if nargin<3 || isempty(Prop)
    Prop={'o','b','b',72,16,2};
end
if nargin<4
    corType='pearson';
end
if nargin<5
    newF=0;
end
% plot scatter
if newF
    figure;
end
format_fig;

% axis(oldaxis);
hold on;
% correlation
if size(X,1)<size(X,2)
    X=X';
end
if size(Y,1)<size(Y,2)
    Y=Y';
end
disc=isnan(X) | isnan(Y);
X(disc)=[];
Y(disc)=[];

bins=prctile(X',Prop{5}:Prop{5}:100);
Ybin(1)=mean(Y(X<=bins(1)));
Xbin(1)=mean(X(X<=bins(1)));
Ybin_sem(1)=sem(Y(X<=bins(1)));
line([1 1]*Xbin(1),[-1 1]*Ybin_sem(1)+Ybin(1),'Color',Prop{2},'LineWidth',Prop{6})
scatter(Xbin(1),Ybin(1),'Marker',Prop{1},'SizeData',Prop{4},'MarkerFaceColor',Prop{3},'MarkerEdgeColor',Prop{2},'LineWidth',Prop{6})
for nbin=2:length(bins)-2
    Ybin(nbin)=mean(Y(X>bins(nbin-1) & X<=bins(nbin)));
    Ybin_sem(nbin)=real(sem(Y(X>bins(nbin-1) & X<=bins(nbin))));
    Xbin(nbin)=mean(X(X>bins(nbin-1) & X<=bins(nbin)));
    Xbin_sem(nbin)=real(sem(X(X>bins(nbin-1) & X<=bins(nbin))));
    
    line([1 1]*Xbin(nbin),[-1 1]*Ybin_sem(nbin)+Ybin(nbin),'Color',Prop{2},'LineWidth',Prop{6})
    line([-1 1]*Xbin_sem(nbin)+Xbin(nbin),[1 1]*Ybin(nbin),'Color',Prop{2},'LineWidth',Prop{6})
    scatter(Xbin(nbin),Ybin(nbin),'Marker',Prop{1},'SizeData',Prop{4},'MarkerFaceColor',Prop{3},'MarkerEdgeColor',Prop{2},'LineWidth',Prop{6})
end
if isempty(nbin)
    nbin=1;
end
nbin=nbin+1;
Ybin(nbin)=mean(Y(X>bins(nbin-1) & X<=Inf));
Xbin(nbin)=mean(X(X>bins(nbin-1) & X<=Inf));
Ybin_sem(nbin)=sem(Y(X>bins(nbin-1) & X<=Inf));
line([1 1]*Xbin(nbin),[-1 1]*Ybin_sem(nbin)+Ybin(nbin),'Color',Prop{2},'LineWidth',Prop{6})
scatter(Xbin(nbin),Ybin(nbin),'Marker',Prop{1},'SizeData',Prop{4},'MarkerFaceColor',Prop{3},'MarkerEdgeColor',Prop{2},'LineWidth',Prop{6})

% scatter(allallPower_binned,allallPC_binned)
if ~isempty(corType)
    [b,stats]=robustfit(X,Y);
    plot((-2*max(abs(X)):(4*max(abs(X)))/10000:2*max(abs(X))),(-2*max(abs(X)):(4*max(abs(X)))/10000:2*max(abs(X)))*b(2)+b(1),'Color',Prop{2},'LineWidth',Prop{6},'LineStyle','--');
    % xlim([-2.1 2.1])
    % ylim([-1 1]*0.15)
    [rhopop, pVpop]=corr(X,Y,'type',corType);
    title({sprintf('r=%1.3f, p=%1.3f',rhopop,pVpop)})
end
if min(Xbin)<0
    xlim([1.2*min(Xbin) 1.2*max(Xbin)])
else
    xlim([0.8*min(Xbin) 1.2*max(Xbin)])
end
ylim([min(-1.2*abs(Ybin_sem)+Ybin) max(1.2*abs(Ybin_sem)+Ybin)])

%%


