clear all
close all

%%% add fieldtrip path
addpath /home/jaco/fieldtrip-read-only/preproc/

%%% create random data
data = randn(32,500,100); % 32 channels, 500 samples, 100 trials
data = cumsum(data,2);

data(32,2:500,:) = data(1,1:499,:);

cfg.chan_sel = 1:32;  % compute for all pairs of channels
cfg.data_sel = 1:500; % compute using all samples
cfg.taus     = [1 2 4 8]; % compute for taus 1 2 4 8
cfg.kernel   = 3; % kernel = 3 (3 samples per symbol)
cfg.sf       = 250;  % sampling frequency
cfg.over_trials       = 1;  % sampling frequency

[sym, count, smi, wsmi] = smi_and_wsmi(data, cfg);


%%% The output variables are structured in cells corresponding 
%%% to each tau.
%%%
%%% Variables:
%%%
%%% sym: the symbolic transformation of the time series (carefull is zero based, symbols form 0 to 5). 
%%% Structure: channel x symbols x trials
%%%
%%% count: the probability (ocurrence rate) of each symbol.
%%% Structure: channel x symbols x trials
%%%
%%% smi: the symbolic mutual information connectivity 
%%% Structure:
%%% channels x channels x trials, with the connectivity between channels
%%%
%%% wsmi: idem to smi but weighted






