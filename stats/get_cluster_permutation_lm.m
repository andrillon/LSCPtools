function [realpos realneg]=get_cluster_permutation_lm(data,group,formulaNames,formula,montecarloalpha,clusteralpha,npermutation,sTime)

% Input
% - data: cell array (1 cell per condition:  channels * times * subjects
% - montecarloalpha
% - clusteralpha
% - npermutation
% - sTime
%% Bootstraping
all_perms=nan(size(data,1),npermutation);
for nperm=1:npermutation
    all_perms(:,nperm)=randperm(size(data,1));
end

%% Stats
fprintf('Calculating LM...');
%     tic;
rdm = data;
fprintf('... sample %4.0f perm %4.0f',0,0)
for nt=1:size(data,2)
    table=array2table([rdm(:,nt) group],'VariableNames',formulaNames);
    mdl = fitlm(table,formula);
    if nt==1
        rd=nan(size(data,2),size(mdl.Coefficients,1)-1);
        pd=nan(npermutation,size(data,2),size(mdl.Coefficients,1)-1);
        rdpv=nan(size(data,2),size(mdl.Coefficients,1)-1);
        pdpv=nan(npermutation,size(data,2),size(mdl.Coefficients,1)-1);
    end
    rd(nt,:) = double(mdl.Coefficients.tStat(2:size(mdl.Coefficients,1)));
    rdpv(nt,:) = double(mdl.Coefficients.pValue(2:size(mdl.Coefficients,1)));
    
    for nperm=1:npermutation
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... sample %4.0f perm %4.0f',nt,nperm)
        pdm = data(all_perms(:,nperm),:);
        table=array2table([pdm(:,nt) group],'VariableNames',formulaNames);
        mdl = fitlm(table,formula);
    
        pd(nperm,nt,:) = double(mdl.Coefficients.tStat(2:size(mdl.Coefficients,1)));
        pdpv(nperm,nt,:) = double(mdl.Coefficients.pValue(2:size(mdl.Coefficients,1)));
    
    end
end
%     toc;
fprintf('\n')

%% cluster statistics
for k=1:size(rd,2)
    fprintf( 'Computing significance\n');
    [realpos{k}, realneg{k}] = findcluster(rd(:,k)', rdpv(:,k)', clusteralpha);
    
    realpos{k}.pmonte = zeros(size(realpos{k}.tclusters));
    realneg{k}.pmonte = zeros(size(realneg{k}.tclusters));
    
    for isim = 1:npermutation
        [simpos, simneg] =  findcluster(squeeze(pd(isim,:,k))', squeeze(pdpv(isim,:,k))', clusteralpha);
        
        maxval = max(simpos.tclusters);
        if ~isempty(maxval)
            realpos{k}.pmonte = realpos{k}.pmonte + (realpos{k}.tclusters < maxval)./npermutation;
        end
        
        minval = min(simneg.tclusters);
        if ~isempty(minval)
            realneg{k}.pmonte = realneg{k}.pmonte + (realneg{k}.tclusters > minval)./npermutation;
        end
    end
    
    
    pmonte = realpos{k}.pmonte;
    goodc = find(pmonte < montecarloalpha);
    contrast = linspace(.5, 1, length(goodc));
    for i = 1:length(goodc)
        ic = goodc(i);
        samples = realpos{k}.clusters == ic;
        cint = [min(sTime(samples)) max(sTime(samples))];
        [peakv,peaki] = max(rd(samples));
        cintsamples = find(samples);
        peakt = sTime(cintsamples(peaki));
        fprintf('\t effect %g: pos | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', k, pmonte(ic), cint, peakt);
        
    end
    
    
    fprintf('\n');
    
    pmonte = realneg{k}.pmonte;
    goodc = find(pmonte < montecarloalpha);
    contrast = linspace(.5, 1, length(goodc));
    for i = 1:length(goodc)
        ic = goodc(i);
        samples = realneg{k}.clusters == ic;
        cint = [min(sTime(samples)) max(sTime(samples))];
        [peakv,peaki] = min(rd(samples));
        cintsamples = find(samples);
        peakt = sTime(cintsamples(peaki));
        fprintf('\t effect %g: neg | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', k, pmonte(ic), cint, peakt);
        
    end
    fprintf('\n')
end

end

%%%%%%%%%%%%
function [ pos, neg ] = findcluster(d, pV, clusteralpha)

    function res = getcluster(d, ok)
        
        clusters = zeros(size(d));
        ntemp = size(d,2);
        res = cell(1,1);
        
        % first find the clusters
        nclusters = 0;
        cluster = 0;
        for it = 1:ntemp
            if ok(1, it)
                if ~cluster
                    nclusters = nclusters + 1;
                    cluster = nclusters;
                end
            else
                cluster = 0;
            end
            
            clusters(1,it) = cluster;
        end
        
        % then compute a sumary, statistics, etc
        iecluster = struct;
        iecluster.clusters = clusters(1,:);
        
        fclusters = zeros(1,nclusters);
        nclusters = max(iecluster.clusters);
        for ic =1:nclusters
            fclusters(ic) = sum(d(1,iecluster.clusters==ic));
        end
        
        iecluster.nclusters = nclusters;
        iecluster.tclusters = fclusters;
        
        res = iecluster;
    end


ok = pV < clusteralpha & d>0;
pos = getcluster(d, ok);

ok = pV < clusteralpha & d<0;
neg = getcluster(d, ok);

%     bigger = repmat(rd,[1 1 1000]) > pd;
end


