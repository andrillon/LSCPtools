function [stat] = compute_Topo_clusterPerm_v2(dat1,dat2,time,chanlabels,fsample,clusteralpha,mc_alpha,nperm,layout,layoutName)

%%% PREPARE DATA
% cfg=[];
cfg.layout = layout;
%     layout = ft_prepare_layout(cfg);       %Layout of the electrodes

% build data for field trip functions (should be D.ftraw, don't works for time-frequency)
data_cond1.label = chanlabels;
data_cond1.time = time;
data_cond1.dimord = 'subj_chan_time';
data_cond1.individual = dat1;
data_cond1.elec = layout.cfg.elec;
data_cond1.cfg = cfg;
data_cond1.fsample = fsample;

data_cond2.label = chanlabels;
data_cond2.time = time;
data_cond2.dimord = 'subj_chan_time';
data_cond2.individual = dat2;
data_cond2.elec = layout.cfg.elec;
data_cond2.cfg = cfg;
data_cond2.fsample = fsample;

% run the clustering/permutation statistics
cfg = [];
cfg.channel          = chanlabels;
cfg.latency          = time;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = clusteralpha;
cfg.clusterstatistic = 'maxsum';
%     cfg.clusterstatistic = 'wcm';
%   cfg.wcm_weight       = 1;

% cfg.minnbchan = minnbchan;
cfg.minnbchan        = 0;
cfg.tail             = 0;
cfg.clustertail      = 0;
%     cfg.alpha            = 0.025;
cfg.alpha            = mc_alpha;
cfg.correcttail      = 'alpha';
cfg.numrandomization = nperm;

cfg_neighb=[];
cfg_neighb.method = 'tri';
if nargin<10
cfg_neighb.layout=layout;
else
cfg_neighb.layout=layoutName;
end
% cfg_neighb.neighbourdist = 0; % for EGI 65, change to 0. For other nets: check distance
neighbours = ft_prepare_neighbours(cfg_neighb, data_cond1);
cfg.neighbours       = neighbours;

%     [stat] = ft_timelockstatistics(cfg, GA_FIC, GA_FC)
cfg.uvar     = 1;
cfg.ivar     = 2;
subj = size(dat1,1);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.statistic        = 'depsamplesT';
cfg.design   = design;

% data_cond2 = data_cond1;
% %         data_cond2.(zparam) = squeeze(permute(data2,[3,1,2]));
% data_cond2.(zparam) =  calc_zvalues(dat,1);
[stat] = ft_timelockstatistics(cfg, data_cond1,data_cond2);