function data_pupil=get_Tobii_cleanpupil(rawdata_pupil,SR,data_time,rawblinks)

data_pupil=rawdata_pupil;
data_pupil(data_pupil<1)=NaN;
blinks.start=rawblinks(:,2);
blinks.end=rawblinks(:,3);
fprintf('... ... cleaning blinks (NaN replacement)\n')
fprintf('%3.0f%%\n',0)
for nbl=1:length(blinks.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
%     [~,startB]=findclosest(data_time,blinks.start(nbl));
%     [~,endB]=findclosest(data_time,blinks.end(nbl));
    startB=blinks.start(nbl);
    endB=blinks.end(nbl);
    %     data_pupil(startB(1):endB(1))=nan;
    if startB-0.2*SR<1 || endB+0.2*SR>length(data_pupil)
        continue;
    end
    temp_interpolate=[mode(data_pupil(startB-0.2*SR:startB-0.1*SR)) ; mode(data_pupil(endB+0.1*SR:endB+0.2*SR))];
    if sum(~isreal(temp_interpolate))~=0
        continue;
    end
    %     [p,S,mu] =polyfit([0 1],temp_interpolate',1);
    p=[];
    p(2)=temp_interpolate(1);
    p(1)=diff(temp_interpolate);
    temp_interpolated=p(1)*(1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1)+p(2);
    data_pupil(startB(1)-0.1*SR:endB(1)+0.1*SR)=temp_interpolated;
%     
%         figure;
%         plot(rawdata_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','r');
%         hold on
%         plot(data_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','k');
%         pause;
end
fprintf('... ... done\n')
