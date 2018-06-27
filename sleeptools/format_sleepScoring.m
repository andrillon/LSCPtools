function [hypnogram_ts arousal_sampled marousal_sampled]=format_sleepScoring(D,filename)

%% Format sleep scoring
load(filename);
hypnogram_ts=repmat(hypnogram,scoringInfo{3}*D.fsample,1);
hypnogram_ts=reshape(hypnogram_ts,1,numel(hypnogram_ts));
hypnogram_ts(D.nsamples+1:end)=[];

%% Add arousals and MAs
arousal_sampled=zeros(1,length(hypnogram_ts));
marousal_sampled=zeros(1,length(hypnogram_ts));

findArousals=round(arousal.*250);
for nA=1:size(findArousals,1)
    arousal_sampled(findArousals(nA,1):findArousals(nA,2))=1;
end
eventsType={events.type};
eventsTimes=[events.time];

findMAbeg=eventsTimes(match_str(eventsType,'MA beg'))';
findMAend=eventsTimes(match_str(eventsType,'MA end'))';
findMAevbeg=eventsTimes(match_str(eventsType,'MAev beg'))';
findMAevend=eventsTimes(match_str(eventsType,'MAev end'))';

findallMAbeg=sort([findMAbeg ; findMAevbeg]);
findallMAend=sort([findMAend ; findMAevend]);

nAc=0;
for nA=1:min([length(findallMAbeg) length(findallMAend)])
    temp=findallMAend-findallMAbeg(nA);
    temp(temp<0)=NaN;
    thisEnd=find(temp==nanmin(temp));
    if findallMAend(thisEnd)-findallMAbeg(nA)<20 || (nA<length(findallMAbeg) && findallMAend(thisEnd)-findallMAbeg(nA+1)<0)
        nAc=nAc+1;
        findMA(nAc,:)=round([findallMAbeg(nA) findallMAend(thisEnd)]*250);
        if findMA(nAc,2)-findMA(nAc,1)>3*250
            arousal_sampled(findMA(nAc,1):findMA(nAc,2))=1;
        else
            marousal_sampled(findMA(nAc,1):findMA(nAc,2))=1;
        end
    end
end

hypnogram_ts(arousal_sampled==1)=0;