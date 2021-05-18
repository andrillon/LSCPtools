function [data_pupil, filt_data_pupil]=get_EyeLink_cleanpupil(rawdata_pupil,SR,data_time,EL_events)

data_pupil=rawdata_pupil;
% data_pupil(data_pupil<100)=NaN;
blinks=EL_events.Blinks;
fprintf('... ... cleaning blinks (linear interpolation)\n')
fprintf('%3.0f%%\n',0)
for nbl=1:length(blinks.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
    startB=find(data_time==blinks.start(nbl));
    endB=find(data_time==blinks.end(nbl));
    %     data_pupil(startB(1):endB(1))=nan;
    if startB-0.2*SR<1 || endB+0.2*SR>length(data_pupil)
        continue;
    end
%     temp_interpolate=[mode(data_pupil(startB-0.2*SR:startB-0.1*SR)) ; mode(data_pupil(endB+0.1*SR:endB+0.2*SR))];
    temp_interpolate=[data_pupil(startB-0.1*SR) ; data_pupil(endB+0.1*SR)];
    if sum(~isreal(temp_interpolate))~=0
        continue;
    end
    %     [p,S,mu] =polyfit([0 1],temp_interpolate',1);
    p=[];
    p(2)=temp_interpolate(1);
    p(1)=diff(temp_interpolate);
    temp_interpolated=p(1)*(1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1)+p(2);
    data_pupil(startB(1)-0.1*SR:endB(1)+0.1*SR)=temp_interpolated;
    
    %     figure;
    %     plot(rawdata_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','r');
    %     hold on
    %     plot(data_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','k');
    %     pause;
end

% filter before NaN removal
fprintf('... ... smoothing (6Hz low-pass)\n')
filt_data_pupil=lowpass(data_pupil, SR, 6, 4);

% add NaNs for missing data
blinks=EL_events.Blinks;
dur_blinks=(blinks.duration)/SR;
closure=blinks;
closure.start(dur_blinks<2)=[];
closure.end(dur_blinks<2)=[];
closure.duration(dur_blinks<2)=[];
fprintf('... ... cleaning eye closure (NaN replacement)\n')
fprintf('%3.0f%%\n',0)
for nbl=1:length(closure.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(closure.start)*100)
    startB=find(data_time==closure.start(nbl));
    endB=find(data_time==closure.end(nbl));
    %     data_pupil(startB(1):endB(1))=nan;
    if startB-0.2*SR<1 || endB+0.2*SR>length(data_pupil) %|| ((endB-startB)/SR)<2
        continue;
    else
        filt_data_pupil(startB-0.1*SR:endB+0.1*SR)=NaN;
        data_pupil(startB-0.1*SR:endB+0.1*SR)=NaN;
    end
end
fprintf('... ... done\n')
