function [pV hplot]=simpleTplot(x,y,newF,colorF,statsF,lineF,transpF,jbFlag,sthFlag,errFlag,lineWidth)

% [pV hplot]=simpleTplot(x,y,newF,colorF,statsF,lineF,transpF,jbFlag,sthFlag)
%
% INPUTS:
% - x: time vector
% - y: data (either vector (averaged-data) or matrix (non-averaged, will plot the mean+SEM) or cell array (will plot mean+SEM for muliple conditions))
% - newF: 1 (plot on a new window) 0 (keep the same window)
% - colorF: color of the plot. if y is vector or matrix 1 color. Or cell
% array of color
% - statsF: [0]: no stat; [1 alpha Ho]: ttest alpha; [2 montecarloalpha
% clusteralpha nperm]: cluster perm ; [3 alpha]: FDR correction
% - lineF: line type of the mean
% - transpF: transparency of the std fill (0->1)
% - jbFlag: use jbfill or not
% - sthFlag: 0 no smooth or smooth with window=sthFlag
% - errFlag: display or not the error estimate (1: yes; 2: dotted line; 0:
% no)
% - lineWidth: width of the lines plotted (normal: 1)
% OUTPUTS:
% - pV: pValues of the t-test
% - hplot: handles of the curves
%

% Thomas Andrillon
% thomas.andrillon@ens.fr

if nargin<3
    newF=1;
end
if nargin<4
    colorF='b';
end
if nargin<5
    statsF=0;
end
if statsF(1)==1 && numel(statsF)==1
    statsF(2)=0.05;
elseif statsF(1)==1 && numel(statsF)==1
    statsF=[2 0.05 0.05 200];
end
if nargin<6
    lineF='-';
end
if nargin<7
    transpF=0.5;
end
if nargin<8
    jbFlag=1;
end
if nargin<9
    sthFlag=0;
end
if nargin<10
    errFlag=1;
end
if nargin<11
    lineWidth=3;
end

if statsF(1)==1 && length(statsF)<3
    statsF(3)=0;
end
pV=[];
hplot=[];
if ~iscell(y)
    if ndims(y)==2
        if newF
            figure;
        end
        if size(y,1)==1
            hold on;
            if sthFlag~=0
                y=fastsmooth(y,sthFlag,3,1);
            end
            plot(x,y,'LineWidth',4,'Color',colorF,'LineStyle',lineF,'LineWidth',lineWidth);
            
        else
            hold on;
            
            y(sum(isnan(y),2)~=0,:)=[];
            yo=y;
            if sthFlag~=0
                for n=1:size(y,1)
                    y(n,:)=fastsmooth(y(n,:),sthFlag,3,1);
                end
            end
            if jbFlag
                jbfill(x,mean(y)+std(y)/sqrt(size(y,1)-1),mean(y)-std(y)/sqrt(size(y,1)-1),colorF,colorF,1,transpF);
                hold on;
                hplot=plot(x,mean(y),'LineWidth',4,'Color',colorF,'LineStyle',lineF);
            else
                if errFlag==1
                    hplot=plot(x,mean(y),'LineWidth',lineWidth+1,'Color',colorF,'LineStyle',lineF); hold on;
                    plot(x,mean(y)+std(y)/sqrt(size(y,1)-1),'LineWidth',lineWidth,'Color',colorF,'LineStyle','--');
                    plot(x,mean(y)-std(y)/sqrt(size(y,1)-1),'LineWidth',lineWidth,'Color',colorF,'LineStyle','--');
                elseif errFlag==2
                    errorbar(x,mean(y),mean(y)+std(y)/sqrt(size(y,1)-1),mean(y)-std(y)/sqrt(size(y,1)-1),'LineWidth',lineWidth,'Color',colorF,'LineStyle',lineF);
                else
                    hplot=plot(x,mean(y),'LineWidth',lineWidth+1,'Color',colorF,'LineStyle',lineF); hold on;
                end
            end
            
            if statsF(1)==1
                [h pV]=ttest(yo,statsF(3),statsF(2),[],1);
                if ~isempty(find(h))
                    myaxis=axis;
                    scatter(x(find(h)),statsF(3)*ones(1,sum(h))+0.1*rand*myaxis(4),'Marker','*','MarkerFaceColor',colorF,'MarkerEdgeColor',colorF,'LineWidth',lineWidth);
                end
            elseif statsF(1)==2
                myaxis=axis;
                tempS{1}(1,:,:)=yo';
                [realpos realneg]=get_cluster_permutation(tempS,statsF(2),statsF(3),statsF(4),x,1);
                pV.realpos=realpos;
                pV.realneg=realneg;
                % plot positive cluster
                for nC=1:realpos{1}.nclusters
                    if realpos{1}.pmonte(nC)<statsF(2)
                        hold on;
                        plot(x(find(realpos{1}.clusters==nC)),1.1*max(mean(y)+std(y)/sqrt(size(y,1)-1))*ones(1,length(find(realpos{1}.clusters==nC))),'Color',colorF,'LineWidth',lineWidth-1,'LineStyle',lineF)
                    end
                end
                % plot negative cluster
                for nC=1:realneg{1}.nclusters
                    if realneg{1}.pmonte(nC)<statsF(2)
                        hold on;
                        plot(x(find(realneg{1}.clusters==nC)),1.1*min(mean(y)-std(y)/sqrt(size(y,1)-1))*ones(1,length(find(realneg{1}.clusters==nC))),'Color',colorF,'LineWidth',lineWidth-1,'LineStyle',lineF)
                    end
                end
            elseif statsF(1)==3
                [h pV]=ttest(yo,statsF(3),statsF(2),[],1);
                [p_fdr, p_masked] = fdr( pV, statsF(2));
                h=zeros(size(pV));
                h(pV<p_fdr)=1;
                if ~isempty(find(h))
                    myaxis=axis;
                    scatter(x(find(h)),zeros(1,sum(h))+1.1*min(mean(y)-std(y)/sqrt(size(y,1)-1)),'Marker','*','MarkerFaceColor',colorF,'MarkerEdgeColor',colorF,'LineWidth',lineWidth);
                end
            end
        end
    else
        error('y needs to be a vector or a 2D matrix!')
    end
elseif iscell(y)
    if newF
        figure;
    end
    for nc=1:length(y)
        ytemp=y{nc};
        if ndims(ytemp)==2
            
            if size(ytemp,1)==1
                hold on;
                if sthFlag~=0
                    ytemp=fastsmooth(ytemp,sthFlag,3,1);
                end
                plot(x,ytemp);
            else
                hold on;
                
                ytemp(sum(isnan(ytemp),2)~=0,:)=[];
                yo=ytemp;
                if sthFlag~=0
                    for n=1:size(ytemp,1)
                        ytemp(n,:)=fastsmooth(ytemp(n,:),sthFlag,3,1);
                    end
                end
                if jbFlag
                    jbfill(x,mean(ytemp)+std(ytemp)/sqrt(size(ytemp,1)-1),mean(ytemp)-std(ytemp)/sqrt(size(ytemp,1)-1),colorF(nc,:),colorF(nc,:),1,transpF);
                    hold on;
                    hplot(end+1)=plot(x,mean(ytemp),'LineWidth',lineWidth-1,'Color',colorF(nc,:),'LineStyle',lineF{nc});
                else
                    hplot(end+1)=plot(x,mean(ytemp),'LineWidth',lineWidth-1,'Color',colorF(nc,:),'LineStyle',lineF{nc}); hold on;
                    plot(x,mean(ytemp)+std(ytemp)/sqrt(size(ytemp,1)-1),'LineWidth',1,'Color',colorF(nc,:),'LineStyle','--');
                    plot(x,mean(ytemp)-std(ytemp)/sqrt(size(ytemp,1)-1),'LineWidth',1,'Color',colorF(nc,:),'LineStyle','--');
                end
            end
        else
            error('y needs to be a vector or a 2D matrix!')
        end
    end
    if statsF(1)~=0
        if length(y)==2
            if statsF(3)==1
            [h pV]=ttest(y{1},y{2}, statsF(2),[],1);
            elseif statsF(3)==2
            [h pV]=ttest2(y{1},y{2}, statsF(2),[],1);
            end
            corrpV=fdr(pV,statsF(2));
            h=zeros(1,length(pV));
            h(pV<corrpV)=1;
            if ~isempty(find(h))
                myaxis=axis;
                scatter(x(find(h)),zeros(1,sum(h))+0.1*rand*myaxis(4),'Marker','*','MarkerFaceColor',colorF(nc,:),'MarkerEdgeColor',colorF(nc,:));
            end
        else
            atemp=[];
            agroup=[];
            for nT=1:length(x)
                for nC=1:length(y)
                    atemp=[atemp ; y{nC}(:,nT)];
                    agroup=[agroup ; nC*ones(length(y{nC}(:,nT)),1)];
                end
                pV(nT) = kruskalwallis(atemp,agroup,'off');
            end
            h=pV<statsF(2);
            if ~isempty(find(h))
                myaxis=axis;
                scatter(x(find(h)),zeros(1,sum(h))+0.1*rand*myaxis(4),'Marker','*','MarkerFaceColor',colorF(nc,:),'MarkerEdgeColor',colorF(nc,:));
            end
        end
    end
end