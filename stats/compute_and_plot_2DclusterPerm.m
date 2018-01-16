function stat=compute_and_plot_2DclusterPerm(data,timeP,freqP,clusteralpha,montecarloalpha,nperm,plotFlag,tailFlag)
% Inputs:
% - data: Rep x time x freq
% - timeP, freqP: time and freq vector of interest

if nargin<8
    tailFlag=0;
end
nSess = size(data, 1);
nTime = size(data, 3);
nFreq = size(data, 2);

% compute t-value
tval = calc_zvalues(data, 1); % ./ (std(data, 1)/sqrt(nSess));
tval(isnan(tval))=0;
% figure; set(gcf, 'Name', 'T-map'), set(gcf, 'Color', [1 1 1])
% imagesc(timeP,freqP,squeeze(tval)')
% set(gca,'Ydir','normal')
% RESHAPE DATA
% ============
dat = zeros(nTime*nFreq, nSess);
for iSess = 1:nSess
    tmp = data(iSess, :, :);
    dat(:, iSess) = squeeze(tmp(:));
end
% the design is data against zero to compute a t-test
dat = [dat, zeros(size(dat))];
design  = [[ones(nSess, 1); 2*ones(nSess, 1)] ,...
    [1:nSess, 1:nSess]'];

% COMPUTE CLUSTER STAT
% ====================
cfg                     = [];
cfg.neighbours          = [];
cfg.dimord              = 'rpt_freq_times';
cfg.freq                = freqP;
cfg.time                = timeP;
cfg.channel             = {'eeg'};
cfg.dim                 = [1 nFreq nTime];

cfg.numrandomization    = nperm;                   % number of randomizations, can be 'all'
cfg.correctm            = 'cluster';              % apply multiple-comparison correction, 'no', 'max', cluster', 'bonferoni', 'holms', 'fdr' (default = 'no')
cfg.alpha               = montecarloalpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.clusteralpha               = montecarloalpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.tail                = tailFlag;                      % -1, 1 or 0 (default = 0)
cfg.correcttail         = 'no';                % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
cfg.ivar                = 1;                      % number or list with indices, independent variable(s)
cfg.uvar                = 2;                      % number or list with indices, unit variable(s)
cfg.wvar                = [];                     % number or list with indices, within-cell variable(s)
cfg.cvar                = [];                     % number or list with indices, control variable(s)
cfg.feedback            = 'textbar';              % 'gui', 'text', 'textbar' or 'no' (default = 'text')
cfg.randomseed          = 'yes';                  % 'yes', 'no' or a number (default = 'yes')

cfg.clusterstatistic    = 'maxsum';               % no, max, maxsize, maxsum, wcm init: maxsum
cfg.clusterthreshold    = 'nonparametric_common'; % parametric, nonparametric_individual, nonparametric_common
cfg.clusteralpha        = clusteralpha;                   %
cfg.clustercritval      = [] ;                    %
cfg.clustertail         = cfg.tail;               %

cfg.statistic           = 'depsamplesT';
cfg.avgoverchan         = 'no';
[stat, cfg]             = ft_statistics_montecarlo(cfg, dat, design);


% PLOT RESULT
% ===========
fprintf('\n****\n')

if plotFlag
    figure; set(gcf, 'Color', [1 1 1]); set(gcf, 'Name', 'Cluster corrected of t-values')
    h = imagesc(timeP,freqP,reshape(stat.stat, [nTime nFreq])');
    newmask = 0.5*double(reshape(stat.mask, [nTime nFreq]))';
    newmask(newmask(:) == 0) = 0;
    alpha(h, newmask)
    title('Significant cluster (opaque)')
    xlabel('Time')
    ylabel('Frequency')
    set(gca,'Ydir','normal')
end

if isfield(stat,'posclusters') && ~isempty([stat.posclusters])
    probClus=[stat.posclusters.prob];
    sigClus=probClus(probClus<montecarloalpha);
    if ~isempty(sigClus)
        for nClus=1:length(sigClus)
            fprintf('... 1 positive cluster, p: %1.3f\n',sigClus(nClus))
        end
    else
        fprintf('... NO significant positive cluster, all p>%g\n',montecarloalpha)
    end
else
    fprintf('... NO positive cluster, all p_alpha>%g\n',clusteralpha)
end

fprintf('****\n')
if isfield(stat,'negclusters') && ~isempty([stat.negclusters])
    probClus=[stat.negclusters.prob];
    sigClus=probClus(probClus<montecarloalpha);
    if ~isempty(sigClus)
        for nClus=1:length(sigClus)
            fprintf('... 1 negative cluster, p: %1.3f\n',sigClus(nClus))
        end
    else
        fprintf('... NO significant negative cluster, all p>%g\n',montecarloalpha)
    end
else
    fprintf('... NO negative cluster, all p_alpha>%g\n',clusteralpha)
end