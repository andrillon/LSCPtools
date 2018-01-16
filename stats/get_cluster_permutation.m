function [realpos realneg]=get_cluster_permutation(data,montecarloalpha,clusteralpha,npermutation,sTime,averageelecFlag,dataPerm)

% Input
% - data: cell array (1 cell per condition:  channels * times * subjects
% - montecarloalpha
% - clusteralpha
% - npermutation
% - sTime
% - averageelecFlag

    Ncond=size(data,2);
if nargin<7
    dataPerm=cell(1,length(data));
end
%% Bootstrap
% fprintf('Calculating permutations...');
for iCond=1:Ncond

    nelec=size(data{iCond},1);
    ntime=size(data{iCond},2);
    nsuj=size(data{iCond},3);
    df=nsuj-1;
    if ~isempty(~dataPerm{iCond})
        bdata{iCond}=dataPerm{iCond};
        npermutation=size(bdata{iCond},4);
    else
        bdata{iCond} = zeros(nelec, ntime, nsuj, npermutation);
        %     tic;
        for ib =1:npermutation
            % permute random number of subjects
            p = 2*(rand(1,nsuj) > .5)-1;
            bp = (permute(repmat(p, [nelec 1 ntime]), [1 3 2]));
            
            bdata{iCond}(:,:,:,ib) = bp .* data{iCond};
        end
    end
    %     toc;
end
fprintf('\n')
%% Stats
rd = cell(1,Ncond);
pd = cell(1,Ncond);
for iCond = 1:Ncond
%     fprintf('Calculating t-values for condition %g...', iCond);
    %     tic;
    rdm = data{iCond};
    pdm = bdata{iCond}(:,:,:,:);
    
    if averageelecFlag
        rdm = mean(rdm,1);
        pdm = mean(pdm,1);
    end
    
    [~,~,~,STATS] = ttest(rdm, 0, clusteralpha, 'both', 3);
    rd{iCond} = STATS.tstat;
    
    [~,~,~,STATS] = ttest(pdm, 0, clusteralpha, 'both', 3);
    pd{iCond} = STATS.tstat;
    %     toc;
end
fprintf('\n')

%% cluster statistics
for iCond = 1:Ncond
        nsuj=size(data{iCond},3);
    df=nsuj-1;
    
%     fprintf( 'Computing significance for %g\n', iCond);
    [realpos, realneg] = findcluster(rd{iCond}, df, clusteralpha);
    
    nelec = size(rd{iCond}, 1);
    
    for ie = 1:nelec
        realpos{ie}.pmonte = zeros(size(realpos{ie}.tclusters));
        realneg{ie}.pmonte = zeros(size(realneg{ie}.tclusters));
    end
    
    for isim = 1:npermutation
        [simpos, simneg] =  findcluster(pd{iCond}(:,:,1,isim), df, clusteralpha);
        for ie = 1:nelec
            
            maxval = max(simpos{ie}.tclusters);
            if ~isempty(maxval)
                realpos{ie}.pmonte = realpos{ie}.pmonte + (realpos{ie}.tclusters < maxval)./npermutation;
            end
            
            minval = min(simneg{ie}.tclusters);
            if ~isempty(minval)
                realneg{ie}.pmonte = realneg{ie}.pmonte + (realneg{ie}.tclusters > minval)./npermutation;
            end
        end
    end
    
    
    for ie = 1:nelec
        pmonte = realpos{ie}.pmonte;
        goodc = find(pmonte < montecarloalpha);
        contrast = linspace(.5, 1, length(goodc));
        for i = 1:length(goodc)
            ic = goodc(i);
            samples = realpos{ie}.clusters == ic;
            cint = [min(sTime(samples)) max(sTime(samples))];
            [peakv,peaki] = max(rd{iCond}(ie,samples));
            cintsamples = find(samples);
            peakt = sTime(cintsamples(peaki));
            fprintf('\tcond %d elec %d | pos | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', iCond,ie, pmonte(ic), cint, peakt);
            
        end
    end
    
%     fprintf('\n');
    
    for ie = 1:nelec
        pmonte = realneg{ie}.pmonte;
        goodc = find(pmonte < montecarloalpha);
        contrast = linspace(.5, 1, length(goodc));
        for i = 1:length(goodc)
            ic = goodc(i);
            samples = realneg{ie}.clusters == ic;
            cint = [min(sTime(samples)) max(sTime(samples))];
            [peakv,peaki] = min(rd{iCond}(ie,samples));
            cintsamples = find(samples);
            peakt = sTime(cintsamples(peaki));
            fprintf('\tcond %d elec %d | neg | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', iCond,ie, pmonte(ic), cint, peakt);
            
        end
    end
    
end
% fprintf('\n')

end

%%%%%%%%%%%%
function [ pos, neg ] = findcluster(d, df, clusteralpha)

    function res = getcluster(d, ok)
        
        clusters = zeros(size(d));
        nelec = size(d,1);
        ntemp = size(d,2);
        res = cell(1,nelec);
        
        % for each electrode...
        for ie = 1:nelec
            
            % first find the clusters
            nclusters = 0;
            cluster = 0;
            for it = 1:ntemp
                if ok(ie, it)
                    if ~cluster
                        nclusters = nclusters + 1;
                        cluster = nclusters;
                    end
                else
                    cluster = 0;
                end
                
                clusters(ie,it) = cluster;
            end
            
            % then compute a sumary, statistics, etc
            iecluster = struct;
            iecluster.clusters = clusters(ie,:);
            
            tclusters = zeros(1,nclusters);
            nclusters = max(iecluster.clusters);
            for ic =1:nclusters
                tclusters(ic) = sum(d(ie,iecluster.clusters==ic));
            end
            
            iecluster.nclusters = nclusters;
            iecluster.tclusters = tclusters;
            
            res{ie} = iecluster;
        end
        
    end


maxt = tinv(1-clusteralpha, df);
ok = d > maxt;
pos = getcluster(d, ok);

mint = tinv(clusteralpha, df);
ok = d < mint;
neg = getcluster(d, ok);

%     bigger = repmat(rd,[1 1 1000]) > pd;
end


