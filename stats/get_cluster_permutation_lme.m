function [realpos realneg]=get_cluster_permutation_lme(data,group,formulaNames,formula,montecarloalpha,clusteralpha,npermutation,sTime)

% Input
% - data: cell array (1 cell per condition:  channels * times * subjects
% - montecarloalpha
% - clusteralpha
% - npermutation
% - sTime

%% Data
ntime=size(table,2);


%% Stats
fprintf('Calculating LME...');
%     tic;
rdm = data;

fprintf('... sample %4.0f perm %4.0f',0,0)
for nt=1:size(data,2)
    table=array2table([rdm(:,nt) group],'VariableNames',formulaNames);
    mdl = fitlme(table,formula);
    if nt==1
        rd=nan(size(data,2),2);
        pd=nan(npermutation,size(data,2),2);
    end
    if strcmp(model,'full')
        rd(nt,:) = cell2mat(TAB(2:4,6));
        rdpv(nt,:) = cell2mat(TAB(2:4,7));
    else
        rd(nt,:) = cell2mat(TAB(2:3,6));
        rdpv(nt,:) = cell2mat(TAB(2:3,7));
    end
    
    for nperm=1:npermutation
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... sample %4.0f perm %4.0f',nt,nperm)
        pdm = data(randperm(size(data,1)),:);
        [~,TAB,~] = anovan(pdm(:,nt), group,'model','full','display','off');
        if strcmp(model,'full')
            pd(nperm,nt,:) = cell2mat(TAB(2:4,6));
            pdpv(nperm,nt,:) = cell2mat(TAB(2:4,7));
        else
            pd(nperm,nt,:) = cell2mat(TAB(2:3,6));
            pdpv(nperm,nt,:) = cell2mat(TAB(2:3,7));
        end
    end
end
%     toc;
fprintf('\n')

%% cluster statistics
for k=1:size(rd,2)
    fprintf( 'Computing significance\n');
    [realpos{k}, realneg{k}] = findcluster(rd(:,k)', df1, df2, clusteralpha);
    
    realpos{k}.pmonte = zeros(size(realpos{k}.tclusters));
    realneg{k}.pmonte = zeros(size(realneg{k}.tclusters));
    
    for isim = 1:npermutation
        [simpos, simneg] =  findcluster(squeeze(pd(isim,:,k))', df1, df2, clusteralpha);
        
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
        fprintf('\t pos | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', pmonte(ic), cint, peakt);
        
    end
    
    
    fprintf('\n');
    
    %         pmonte = realneg{k}.pmonte;
    %         goodc = find(pmonte < montecarloalpha);
    %         contrast = linspace(.5, 1, length(goodc));
    %         for i = 1:length(goodc)
    %             ic = goodc(i);
    %             samples = realneg{k}.clusters == ic;
    %             cint = [min(sTime(samples)) max(sTime(samples))];
    %             [peakv,peaki] = min(rd(samples));
    %             cintsamples = find(samples);
    %             peakt = sTime(cintsamples(peaki));
    %             fprintf('\t neg | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', pmonte(ic), cint, peakt);
    %
    %         end
    %         fprintf('\n')
end

end

%%%%%%%%%%%%
function [ pos, neg ] = findcluster(d, df1, df2, clusteralpha)

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


maxt = finv(1-clusteralpha, df1, df2);
ok = d > maxt;
pos = getcluster(d, ok);

mint = finv(clusteralpha, df1, df2);
ok = d < mint;
neg = getcluster(d, ok);

%     bigger = repmat(rd,[1 1 1000]) > pd;
end


