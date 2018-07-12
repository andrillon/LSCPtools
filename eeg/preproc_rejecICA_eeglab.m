function Dnew=preproc_rejecICA_eeglab(D,param)

ALLEEG = eeglab;

spm_datapath = param.data_path;
file_name = D.fname;
ica_file_name = [file_name(1:end-4) '_ica.set'];
fprintf('\n\nProcessing %s...\n',file_name);

% load ICA 
EEG = pop_loadset( ica_file_name, spm_datapath);

% update channel location (just in case)
netfile = [param.path_locfile]; % for wanderlust .xyz
EEG = pop_chanedit(EEG,'load',{netfile, 'filetype', 'autodetect'});

% reject channels that were not included in the ICA composition
ch_usedInICA=EEG.icachansind;
EEG = pop_select( EEG,'nochannel',setdiff(1:size(EEG.data,1),EEG.icachansind));


% Visually inspect the data, and select components for removal and trials to reject
compbad=[];
if strcmp(param.ica.inspecmethod,'adjust')
    [art, horiz, vert, blink, disc, ...
        soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin] = ADJUST (EEG,['adjust_report_' file_name(1:end-4) '.txt']);
    
    % TODO The best would be to do it afterwards and contrast
    % automatic and manually select components, but
    % pop_selectcomps_ADJ is not returning the automatic ones
    % in the EEG object (the highlited in red), so I put all of
    % them as selected so we can unselect them if needed
    EEG.reject.gcompreject(art) = 1;
    compbad = find(EEG.reject.gcompreject);
    
elseif strcmp(param.ica.inspecmethod,'visual')
    EEG = pop_selectcomps(EEG,1:size(EEG.icaweights,1));
    f = warndlg('Click here when finished COMPONENT selection.');
            waitfor(f);
    compbad = find(EEG.reject.gcompreject);
end

if ~isempty(compbad)
    fprintf('... %g components selected for rejection via %s\n',length(compbad),param.ica.inspecmethod)
    save([spm_datapath filesep 'rejcomp_ica_' param.ica.inspecmethod '_'  file_name(1:end-4) ],'compbad')
    
    EEG = pop_subcomp(EEG, compbad);
    
    S = [];
    S.D = D;
    new_smp_fname = fullfile(spm_datapath,['ICA_' param.ica.inspecmethod '_' D.fname]);
    S.outfile = new_smp_fname;
    Dnew = spm_eeg_copy(S);
    
    % updating meeg-object with new data (finally...)
    Dnew(ch_usedInICA,:,:) = EEG.data; 
    Dnew.save;
end

