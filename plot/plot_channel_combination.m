% Create paired electrodes variable
% Check if appariment is good

function plot_channel_combination(net_name,path_name,ElecToHighlight)

if isempty(net_name); net_name='egi256_GSN_HydroCel'; end
if isempty(path_name); path_name='/media/Data/Work/Scripts/newSPM/projects/SleepLexicalDecision'; end
if nargin<3, ElecToHighlight=[]; end
if isnumeric(ElecToHighlight) && strcmp(net_name(1:3),'egi')
    NumElecToHighlight=ElecToHighlight;
    ElecToHighlight=[];
    for n=1:length(NumElecToHighlight)
        ElecToHighlight{n}=sprintf('E%g',NumElecToHighlight(n));
    end
end

    %% Load coordinates
    fid = fopen(fullfile(path_name,[net_name '.sfp']),'r');
    sensC = textscan(fid,'%s %f %f %f',[4 inf]);
    fclose(fid);

    
    nchannels = size(sensC{1}(4:end),1);
    sens.pnt = [sensC{2}(4:end) sensC{3}(4:end) sensC{4}(4:end)];
    sens.label = sensC{1}(4:end);
    for n=1:nchannels
        sens.label{n}(1)='E';
    end
    sens.unit='mm';
    [xy,label] = spm_eeg_project3D(sens, 'EEG'); sens.pnt=xy';
    
    h=figure;
        set(h,'Position',[427   427   250   235],'Color','w')
    plot(sens.pnt(:,1),sens.pnt(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k')
    
    hold on;
    for n=1:nchannels;
        if ismember(sens.label{n},ElecToHighlight)
            plot(sens.pnt(n,1),sens.pnt(n,2),'o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r')
        end
    end
    set(gca,'box','off','XColor','w','YColor','w')
    
