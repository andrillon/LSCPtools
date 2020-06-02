function D=preproc_runICA_eeglab(param)
% Original from Leonardo Barbosa

if ~isfield(param.ica, 'interactive'), param.ica.interactive = false; end
if ~isfield(param.ica, 'useadjust'), param.ica.useadjust = false; end



% select the correct file
spm_datapath = param.data_path;
file_name = param.file_name;
fprintf('\n\nProcessing %s...\n',file_name);

% load spm data
D = spm_eeg_load([spm_datapath filesep file_name]);
%     rmpath(genpath(param.spm12_path));

% restart eeglab
ALLEEG = eeglab;

% Run the ICA (or load if dataset already exists)
ica_file_name = [file_name(1:end-4) '_ica.set'];

% convert SPM data to EEGLAB
EEG = pop_fileio([spm_datapath filesep file_name]);
% EEG = pop_select( EEG,'nochannel',setdiff(1:D.nchannels,match_str(D.chantype,'EEG')));



netfile = [param.name_sensfile]; % for wanderlust .xyz
% if strcmp(param.type_sensfile,'.xyz')
%     EEG = pop_chanedit(EEG,'load',{netfile, 'filetype', param.type_sensfile});
% else
EEG = pop_chanedit(EEG,'load',{netfile, 'filetype', 'autodetect'});
% end
% EEG = pop_select( EEG,'nochannel',badChannels);

EEG = eeg_checkset( EEG );
eeglab redraw
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, 1);

% % Reject bad channels
% fprintf('Computing kurtosis for channels...\n');
% [ measure indelec ] = rejkurt( reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3)), 5, [], 2);
% % reject aberrant channels
% maxVal=squeeze(max(abs(EEG.data),[],2));
% wrongValues=maxVal>300;
% pprtionBad=mean(wrongValues,2);
% badChannels=find(pprtionBad>1/3);
% badChannels=[find(indelec); badChannels];
% fprintf('... found %g bad channels\n',length(badChannels))

% Run the ICA
chanica = match_str(D.chantype,'EEG')';
if strcmp(param.ica.icatype, 'pca')
    EEG = pop_runica(EEG, 'icatype', 'runica', 'pca', param.ica.pcanumcomp, 'chanind', chanica);
elseif strcmp(param.ica.icatype, 'runicanoext')
    EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', chanica);
elseif strcmp(param.ica.icatype, 'runica') % recommended
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'chanind', chanica);
else
    EEG = pop_runica(EEG, 'icatype', param.ica.icatype, 'extended', 1, 'chanind', chanica);
end
EEG = eeg_checkset(EEG);

EEG.setname = ica_file_name;
[ALLEEG,EEG] = eeg_store(ALLEEG,EEG,1);
EEG = pop_saveset(EEG, ica_file_name, spm_datapath);
[ALLEEG,EEG] = eeg_store(ALLEEG, EEG, 1);


%     addpath(genpath(param.spm12_path));
