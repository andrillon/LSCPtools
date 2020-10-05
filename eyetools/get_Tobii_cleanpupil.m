function data_pupil=get_Tobii_cleanpupil(rawdata_pupil,SR,data_time,rawblinks)

data_pupil=rawdata_pupil;
% data_pupil(data_pupil<100)=NaN;
% blinks=EL_events.Blinks;
blinks.start=rawblinks(:,2);
blinks.end=rawblinks(:,3);
blinks.duration=blinks.end-blinks.start;
fprintf('... ... cleaning blinks (linear interpolation)\n')
fprintf('%3.0f%%\n',0)
nbl=0;
while nbl<length(blinks.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
    nbl=nbl+1;
    startB=blinks.start(nbl);
    endB=blinks.end(nbl);
    %     data_pupil(startB(1):endB(1))=nan;
    if startB-0.1*SR<1 || endB+0.1*SR>length(data_pupil)
        continue;
    end
    before_seg=data_pupil(startB-0.1*SR:startB);
    before_seg(before_seg<2)=[];
    after_seg=data_pupil(endB:endB+0.1*SR);
    after_seg(after_seg<2)=[];
    c=1;
     while isempty(before_seg) && nbl-c>2
        c=c+1;
        startB=blinks.start(nbl-c);
        before_seg=data_pupil(startB-0.1*SR:startB);
        before_seg(before_seg<2)=[];
        after_seg=data_pupil(endB:endB+0.1*SR);
        after_seg(after_seg<2)=[];
%         if isempty(before_seg) || isempty(after_seg)
%             pause;
%         end
    end
    while isempty(after_seg) && nbl<length(blinks.start)-1
        nbl=nbl+1;
        endB=blinks.end(nbl);
        before_seg=data_pupil(startB-0.1*SR:startB);
        before_seg(before_seg<2)=[];
        after_seg=data_pupil(endB:min(endB+0.1*SR,length(data_pupil)));
        after_seg(after_seg<2)=[];
%         if isempty(before_seg) || isempty(after_seg)
%             pause;
%         end
    end
    temp_interpolate=[median(before_seg) ; median(after_seg)];
%     if min(temp_interpolate)<2 || sum(isnan(temp_interpolate))~=0
%         pause;
%     end
    if sum(~isreal(temp_interpolate))~=0
        continue;
    end
    
    %     [p,S,mu] =polyfit([0 1],temp_interpolate',1);
    p=[];
    p(2)=temp_interpolate(1);
    p(1)=diff(temp_interpolate);
    temp_interpolated=p(1)*(1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1)+p(2);
    data_pupil(startB(1)-0.1*SR:endB(1)+0.1*SR)=temp_interpolated;
%     if min(temp_interpolated)<2
%         continue;
%     end
    %     figure;
    %     plot(rawdata_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','r');
    %     hold on
    %     plot(data_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','k');
    %     pause;
end

% % filter before NaN removal
% fprintf('... ... smoothing (6Hz low-pass)\n')
% filt_data_pupil=lowpass(data_pupil, SR, 6, 4);

% add NaNs for missing data
% dur_blinks=(blinks.duration)/SR;
% closure=blinks;
% closure.start(dur_blinks<2)=[];
% closure.end(dur_blinks<2)=[];
% closure.duration(dur_blinks<2)=[];
% fprintf('... ... cleaning eye closure (NaN replacement)\n')
% fprintf('%3.0f%%\n',0)
% for nbl=1:length(closure.start)
%     fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(closure.start)*100)
%     startB=closure.start(nbl);
%     endB=closure.end(nbl);
%     %     data_pupil(startB(1):endB(1))=nan;
%     if startB-0.2*SR<1 || endB+0.2*SR>length(data_pupil) %|| ((endB-startB)/SR)<2
%         continue;
%     else
%         filt_data_pupil(startB-0.1*SR:endB+0.1*SR)=NaN;
%         data_pupil(startB-0.1*SR:endB+0.1*SR)=NaN;
%     end
% end
fprintf('... ... done\n')

%
%
% data_pupil=rawdata_pupil;
% data_pupil(data_pupil<1)=NaN;
% blinks.start=rawblinks(:,2);
% blinks.end=rawblinks(:,3);
% fprintf('... ... cleaning blinks (NaN replacement)\n')
% fprintf('%3.0f%%\n',0)
% for nbl=1:length(blinks.start)
%     fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
% %     [~,startB]=findclosest(data_time,blinks.start(nbl));
% %     [~,endB]=findclosest(data_time,blinks.end(nbl));
%     startB=blinks.start(nbl);
%     endB=blinks.end(nbl);
%     %     data_pupil(startB(1):endB(1))=nan;
%     if startB-0.2*SR<1 || endB+0.2*SR>length(data_pupil)
%         continue;
%     end
%     temp_interpolate=[mode(data_pupil(startB-0.2*SR:startB-0.1*SR)) ; mode(data_pupil(endB+0.1*SR:endB+0.2*SR))];
%     if sum(~isreal(temp_interpolate))~=0
%         continue;
%     end
%     %     [p,S,mu] =polyfit([0 1],temp_interpolate',1);
%     p=[];
%     p(2)=temp_interpolate(1);
%     p(1)=diff(temp_interpolate);
%     temp_interpolated=p(1)*(1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1)+p(2);
%     data_pupil(startB(1)-0.1*SR:endB(1)+0.1*SR)=temp_interpolated;
% %
% %         figure;
% %         plot(rawdata_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','r');
% %         hold on
% %         plot(data_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','k');
% %         pause;
% end
% fprintf('... ... done\n')
