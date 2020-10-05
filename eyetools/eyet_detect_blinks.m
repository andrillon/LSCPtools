%%%%%%%%
%%%% Blink detection using the 'validity' measure from the Tobii system
%%%% Will be detected as a blink, episods with eye-closed of a certain
%%%% durations
%%%%
%%%% Inputs:
%%%% - data (2D or 3D matrix with validity of eye-recordings for both eyes
%%%%   for continous (samples*eye) data or epoched (trial*samples*eye))
%%%% - param.mindur (min duration of a blink in seconds)
%%%% - param.maxdur (max duration of a blink in seconds)
%%%% - param.fs     (sampling rate)
%%%% - param.paired (require a blink on both eyes (1) or not (0))
%%%%
%%%% Outputs:
%%%% - blinks (matrix with as many rows as blinks, columns are: (1) beginning
%%%%   (2) end sample (3) flag of detection on eye (0 both, 1 righ-eye, 2 left)

function blinks=eyet_detect_blinks(data,param)
blinks=[];

% if epoched data
if ndims(data)==3
    for i=1:size(data,1)
        temp=squeeze(data(i,:,:));
        if param.paired
            eyeclosed=find(temp(:,1)==4 & temp(:,2)==4)'; % for tobii system, 4 means no data so eye-closed
            eyeclosed=[-Inf eyeclosed Inf];
            findboundaries=find(eyeclosed(2:end)-eyeclosed(1:end-1)~=1);
            epochs=[];
            for j=1:length(findboundaries)-1
                epochs=[epochs ; [i eyeclosed(findboundaries(j)+1) eyeclosed(findboundaries(j+1)) findboundaries(j+1)-findboundaries(j) (findboundaries(j+1)-findboundaries(j))/param.fs 0]];
            end
            if ~isempty(epochs)
                blinks=[blinks ; epochs(epochs(:,5)>param.mindur & epochs(:,5)<param.maxdur,:)];
            end
        else
            for neye=1:2
                eyeclosed=find(temp(:,neye)==4)'; % for tobii system, 4 means no data so eye-closed
                eyeclosed=[-Inf eyeclosed Inf];
                findboundaries=find(eyeclosed(2:end)-eyeclosed(1:end-1)~=1);
                epochs=[];
                for j=1:length(findboundaries)-1
                    epochs=[epochs ; [i eyeclosed(findboundaries(j)+1) eyeclosed(findboundaries(j+1)) findboundaries(j+1)-findboundaries(j) (findboundaries(j+1)-findboundaries(j))/param.fs neye]];
                end
                
                if ~isempty(epochs)
                    blinks=[blinks ; epochs(epochs(:,5)>param.mindur & epochs(:,5)<param.maxdur,:)];
                end
            end
        end
    end
    % if continuous data
else
    temp=squeeze(data(:,:));
    if param.paired
        eyeclosed=find(temp(:,1)==4 & temp(:,2)==4)'; % for tobii system, 4 means no data so eye-closed
        eyeclosed=[-Inf eyeclosed Inf];
        findboundaries=find(eyeclosed(2:end)-eyeclosed(1:end-1)~=1);
        epochs=[];
        for j=1:length(findboundaries)-1
            epochs=[epochs ; [1 eyeclosed(findboundaries(j)+1) eyeclosed(findboundaries(j+1)) findboundaries(j+1)-findboundaries(j) (findboundaries(j+1)-findboundaries(j))/param.fs 0]];
        end
        
        if ~isempty(epochs)
            blinks=[blinks ; epochs(epochs(:,5)>param.mindur & epochs(:,5)<param.maxdur,:)];
        end
    else
        for neye=1:2
            eyeclosed=find(temp(:,neye)==4)'; % for tobii system, 4 means no data so eye-closed
            eyeclosed=[-Inf eyeclosed Inf];
            findboundaries=find(eyeclosed(2:end)-eyeclosed(1:end-1)~=1);
            epochs=[];
            for j=1:length(findboundaries)-1
                epochs=[epochs ; [1 eyeclosed(findboundaries(j)+1) eyeclosed(findboundaries(j+1)) findboundaries(j+1)-findboundaries(j) (findboundaries(j+1)-findboundaries(j))/param.fs neye]];
            end
            
            if ~isempty(epochs)
                blinks=[blinks ; epochs(epochs(:,5)>param.mindur & epochs(:,5)<param.maxdur,:)];
            end
        end
    end
    old_blinks=blinks;
    blinks=[];
    if param.paired==0
        for neye=1:2
            temp_blinks=old_blinks(old_blinks(:,end)==neye,:);
            dist_blinks=temp_blinks(2:end,2)-temp_blinks(1:end-1,3);
            c=1;
            co=0;
            while c<length(dist_blinks)
                if dist_blinks(c)<param.mindist*param.fs
                    temp_blinks(c,:)=[temp_blinks(c,1:2) temp_blinks(c+1,3) temp_blinks(c+1,3)-temp_blinks(c,2) (temp_blinks(c+1,3)-temp_blinks(c,2))/param.fs temp_blinks(c,end)];
                    temp_blinks(c+1,:)=[];
                    dist_blinks=temp_blinks(2:end,2)-temp_blinks(1:end-1,3);
                    co=co+1;
                else
                    c=c+1;
                end
            end
            blinks=[blinks ; temp_blinks];
        end
    end
end