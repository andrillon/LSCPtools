function [saccades_X, saccades_Y]=eyet_detect_saccaces(eyeposX,eyeposY,param)

% set parameters
% param=[];
% param.dist_exclusion=0.02; % 20ms
% param.slope=3; % minimal slope
% param.amplitude=0.05; % minimal slope
% param.fs

% on x-coordinate
segments_X=[];
for j=1:size(eyeposX,2)
    temp=eyeposX(:,j);
    
    segmentsbounds=find(diff(sign(diff(temp))));
    temp_segments=[[1 ; segmentsbounds(1:end-1)+1] [segmentsbounds(1:end)+1]];
    temp_segments=[temp_segments temp(temp_segments(:,1)) temp(temp_segments(:,2)) abs(temp(temp_segments(:,1))-temp(temp_segments(:,2))) abs(temp(temp_segments(:,1))-temp(temp_segments(:,2)))./((temp_segments(:,2)-temp_segments(:,1))/param.fs) j*ones(size(temp_segments,1),1)];
    segments_X=[segments_X ; temp_segments];
    
    %     dir=temp(2)-temp(1);
    %     n=2;
    %     startseg=1;
    %     sizeseg=1;
    %     while n<length(temp)
    %         if sign(temp(n+1)-temp(n))==sign(dir)
    %             n=n+1;
    %             sizeseg=sizeseg+1;
    %         else
    %             endseg=n;
    %             segments_X=[segments_X ; [startseg endseg sizeseg temp(startseg) temp(endseg) abs(temp(startseg)-temp(endseg)) abs(temp(startseg)-temp(endseg))/(sizeseg/param.fs) n]];
    %
    %             dir=temp(n+1)-temp(n);
    %             startseg=n;
    %             sizeseg=1;
    %             n=n+1;
    %         end
    %     end
end

% clean x coordinate
saccades_X=segments_X;
saccades_X(saccades_X(:,end-2)<param.amplitude & saccades_X(:,end-1)<param.slope,:)=[];
saccades_X(find(diff(saccades_X(:,1))/param.fs<param.dist_exclusion),:)=[];
% eliminate if signa lost
flag=[];
for nsac=1:size(saccades_X,1)
    if min(saccades_X(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)))<1 || max(saccades_X(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)))>size(eyeposX,1)
        flag(nsac)=1;
    else
        tp=eyeposX(saccades_X(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)),segments_X(nsac,end));
        flag(nsac)=~isempty(isnan(tp));
    end
end
segments_X(flag,:)=[];

if size(eyeposX,2)==2 && sum(saccades_X(:,end)==1)~=0 && sum(saccades_X(:,end)==2)~=0
    flag=zeros(1,size(saccades_X,1));
    temp_sac=saccades_X(saccades_X(:,end)==1,:);
    temp_sac2=saccades_X(saccades_X(:,end)==2,:);
    idx_sac1=find(saccades_X(:,end)==1);
    idx_sac2=find(saccades_X(:,end)==2);
    for nsac=1:size(temp_sac,1)
        [time2,idx]=findclosest(temp_sac2(:,1),temp_sac(nsac,1));
        if abs(time2-temp_sac(nsac,1))/param.fs>param.dist_exclusion
            flag(idx_sac1(nsac))=1;
        else
            flag(idx_sac1(nsac))=0;
        end
    end
    
    for nsac=1:size(temp_sac2,1)
        [time2,idx]=findclosest(temp_sac(:,1),temp_sac2(nsac,1));
        if abs(time2-temp_sac2(nsac,1))/param.fs>param.dist_exclusion
            flag(idx_sac2(nsac))=1;
        else
            flag(idx_sac2(nsac))=0;
        end
    end
    segments_X(find(flag),:)=[];
elseif size(eyeposX,2)==2 && (sum(segments_X(:,end)==1)==0 || sum(segments_X(:,end)==2)==0)
    segments_X=[];
end

% on Y-coordinate
segments_Y=[];
for j=1:size(eyeposX,2)
    temp=eyeposY(:,j);
    
    segmentsbounds=find(diff(sign(diff(temp))));
    temp_segments=[[1 ; segmentsbounds(1:end-1)+1] [segmentsbounds(1:end)+1]];
    temp_segments=[temp_segments temp(temp_segments(:,1)) temp(temp_segments(:,2)) abs(temp(temp_segments(:,1))-temp(temp_segments(:,2))) abs(temp(temp_segments(:,1))-temp(temp_segments(:,2)))./((temp_segments(:,2)-temp_segments(:,1))/param.fs) j*ones(size(temp_segments,1),1)];
    segments_Y=[segments_Y ; temp_segments];
    
end

% clean x coordinate
saccades_Y=segments_Y;
saccades_Y(saccades_Y(:,end-2)<param.amplitude & saccades_Y(:,end-1)<param.slope,:)=[];
saccades_Y(find(diff(saccades_Y(:,1))/param.fs<param.dist_exclusion),:)=[];
% eliminate if signa lost
flag=[];
for nsac=1:size(saccades_Y,1)
    if min(saccades_Y(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)))<1 || max(saccades_Y(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)))>size(eyeposY,1)
        flag(nsac)=1;
    else
        tp=eyeposY(saccades_Y(nsac,1)+(-round(param.dist_exclusion*param.fs):round(param.dist_exclusion*param.fs)),saccades_Y(nsac,end));
        flag(nsac)=~isempty(isnan(tp));
    end
end
saccades_Y(flag,:)=[];

if size(eyeposY,2)==2 && sum(saccades_Y(:,end)==1)~=0 && sum(saccades_Y(:,end)==2)~=0
    flag=zeros(1,size(saccades_Y,1));
    temp_sac=saccades_Y(saccades_Y(:,end)==1,:);
    temp_sac2=saccades_Y(saccades_Y(:,end)==2,:);
    idx_sac1=find(saccades_Y(:,end)==1);
    idx_sac2=find(saccades_Y(:,end)==2);
    for nsac=1:size(temp_sac,1)
        [time2,idx]=findclosest(temp_sac2(:,1),temp_sac(nsac,1));
        if abs(time2-temp_sac(nsac,1))/param.fs>param.dist_exclusion
            flag(idx_sac1(nsac))=1;
        else
            flag(idx_sac1(nsac))=0;
        end
    end
    
    for nsac=1:size(temp_sac2,1)
        [time2,idx]=findclosest(temp_sac(:,1),temp_sac2(nsac,1));
        if abs(time2-temp_sac2(nsac,1))/param.fs>param.dist_exclusion
            flag(idx_sac2(nsac))=1;
        else
            flag(idx_sac2(nsac))=0;
        end
    end
    
    saccades_Y(find(flag),:)=[];
elseif size(eyeposY,2)==2 && (sum(saccades_Y(:,end)==1)==0 || sum(saccades_Y(:,end)==2)==0)
    saccades_Y=[];
end