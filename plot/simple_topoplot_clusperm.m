function [stat] = simple_topoplot_clusperm(dat1,dat2,time,chanlabels,fsample,clusteralpha,mc_alpha,nperm,layout,layoutName,paths)

%% Do cluster permutation
%%% PREPARE DATA
% cfg=[];
addpath(paths.spm)
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
cfg_neighb.layout=layoutName;
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

%% Summarise clusters
% Positive cluster
if ~isempty(stat.posclusters)
    poscluster_idx=find([stat.posclusters.prob]<mc_alpha);
else
    poscluster_idx=[];
end
if ~isempty(stat.negclusters)
    negcluster_idx=find([stat.negclusters.prob]<mc_alpha);
else
    negcluster_idx=[];
end
fprintf('\nSummary:\n')
fprintf('... found %g positive clusters (p_mc<%g)\n',length(poscluster_idx),mc_alpha);
chan_clustersID=zeros(1,length(chanlabels));
for ncl=1:length(poscluster_idx)
    fprintf('... ... cluster %g: p_mc=%1.4f n=%g elec\n',ncl,stat.posclusters(ncl).prob,sum(stat.posclusterslabelmat==ncl));
    chan_clustersID(find(stat.posclusterslabelmat==ncl))=1;
end
fprintf('... found %g negative clusters (p_mc<%g)\n',length(negcluster_idx),mc_alpha);
for ncl=1:length(negcluster_idx)
    fprintf('... ... cluster %g: p_mc=%1.4f n=%g elec\n',ncl,stat.negclusters(ncl).prob,sum(stat.negclusterslabelmat==ncl));
    chan_clustersID(find(stat.negclusterslabelmat==ncl))=-1;
end

%% Plot topography and cluster
addpath([paths.eeglab filesep 'functions' filesep 'sigprocfunc'])
addpath([paths.eeglab filesep 'functions' filesep 'guifunc'])
addpath([paths.eeglab filesep 'functions' filesep 'adminfunc'])

figure;
temp_toplot=nanmean(dat1-dat2,1);
topoplot(temp_toplot, layout.chaninfo,'style','map','emarker2',{find(chan_clustersID==1 | chan_clustersID==-1 ),'o','k',10,2});
topoplot(chan_clustersID, layout.chaninfo,'style','contour','numcontour',3);
caxis([-1 1]*max(abs(temp_toplot)));
format_fig;

%% Clean paths
rmpath(paths.spm)
rmpath([paths.eeglab filesep 'functions' filesep 'sigprocfunc'])
rmpath([paths.eeglab filesep 'functions' filesep 'guifunc'])
rmpath([paths.eeglab filesep 'functions' filesep 'adminfunc'])