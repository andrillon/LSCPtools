function [realpos realneg]=get_cluster_permutation_aov(data,group,montecarloalpha,clusteralpha,npermutation,sTime,model,contVarIndex)

% Input 
% - data: subjects/conditions * times)
% - group: describes lines
% - montecarloalpha
% - clusteralpha
% - npermutation
% - sTime
% - averageelecFlag

if nargin<7
    model='full';
end
if nargin<7
    contVarIndex=[];
end
%% Bootstrap
% fprintf('Calculating permutations...');
if size(group,2)==1
    if ~iscategorical(group)
        data=data(~isnan(group),:);
        group=group(~isnan(group));
    else
        data=data(~isundefined(group),:);
        group=group(~isundefined(group));
    end
    ntime=size(data,2);
    nsuj=size(data,1);
    df1=length(unique(group))-1;
    df2=nsuj-1;
    
    all_perms=nan(nsuj,npermutation);
    for nperm=1:npermutation
        all_perms(:,nperm)=randperm(size(data,1));
    end

    %
    % bdata = zeros(nsuj, ntime, npermutation);
    % %     tic;
    % for ib =1:npermutation
    %     % permute random number of subjects
    %     p = 2*(rand(1,nsuj) > .5)-1;
    %     bp = (permute(repmat(p, [nsuj 1 ntime]), [1 3 2]));
    %
    %     bdata(:,:,:,ib) = bp .* data;
    % end
    % %     toc;
    % fprintf('\n')
    
    %% Stats
    fprintf('Calculating F-values (anova)...');
    %     tic;
    rdm = data;
    
    for nt=1:size(data,2)
        [~,TAB,~] = anova1(rdm(:,nt), group,'off');
        rd(nt) = TAB{2,5};
        
        for nperm=1:npermutation
            pdm = data(all_perms(:,nperm),:);
            [~,TAB,~] = anova1(pdm(:,nt), group,'off');
            pd(nperm,nt) = TAB{2,5};
        end
    end
    %     toc;
    fprintf('\n')
    
    %% cluster statistics
    fprintf( 'Computing significance\n');
    [realpos, realneg] = findcluster(rd, df1, df2, clusteralpha);
    
    realpos.pmonte = zeros(size(realpos.tclusters));
    realneg.pmonte = zeros(size(realneg.tclusters));
    
    for isim = 1:npermutation
        [simpos, simneg] =  findcluster(pd(isim,:), df1, df2, clusteralpha);
        
        maxval = max(simpos.tclusters);
        if ~isempty(maxval)
            realpos.pmonte = realpos.pmonte + (realpos.tclusters < maxval)./npermutation;
        end
        
        minval = min(simneg.tclusters);
        if ~isempty(minval)
            realneg.pmonte = realneg.pmonte + (realneg.tclusters > minval)./npermutation;
        end
    end
    
    
    pmonte = realpos.pmonte;
    goodc = find(pmonte < montecarloalpha);
    contrast = linspace(.5, 1, length(goodc));
    for i = 1:length(goodc)
        ic = goodc(i);
        samples = realpos.clusters == ic;
        cint = [min(sTime(samples)) max(sTime(samples))];
        [peakv,peaki] = max(rd(samples));
        cintsamples = find(samples);
        peakt = sTime(cintsamples(peaki));
        fprintf('\t pos | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', pmonte(ic), cint, peakt);
        
    end
    
    
%     fprintf('\n');
%     
%     pmonte = realneg.pmonte;
%     goodc = find(pmonte < montecarloalpha);
%     contrast = linspace(.5, 1, length(goodc));
%     for i = 1:length(goodc)
%         ic = goodc(i);
%         samples = realneg.clusters == ic;
%         cint = [min(sTime(samples)) max(sTime(samples))];
%         [peakv,peaki] = min(rd(samples));
%         cintsamples = find(samples);
%         peakt = sTime(cintsamples(peaki));
%         fprintf('\t neg | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', pmonte(ic), cint, peakt);
%         
%     end
    
    fprintf('\n')
    
else
    
    ntime=size(data,2);
    nsuj=size(data,1);
    df1=length(unique(group))-1;
    df2=nsuj-1;
    all_perms=nan(nsuj,npermutation);
    for nperm=1:npermutation
        all_perms(:,nperm)=randperm(size(data,1));
    end
    
    %% Stats
    fprintf('Calculating F-values (anova)...');
    %     tic;
    rdm = data;
    
    if strcmp(model,'full')
        rd=nan(size(data,2),3);
        pd=nan(npermutation,size(data,2),3);
    else
        rd=nan(size(data,2),2);
        pd=nan(npermutation,size(data,2),2);
    end
    fprintf('... sample %4.0f perm %4.0f',0,0)
    for nt=1:size(data,2)
             [~,TAB,~] = anovan(rdm(:,nt), group,'model',model,'display','off','continuous',contVarIndex);
       if strcmp(model,'full')
            rd(nt,:) = cell2mat(TAB(2:4,6));
            rdpv(nt,:) = cell2mat(TAB(2:4,7));
        else
            rd(nt,:) = cell2mat(TAB(2:3,6));
            rdpv(nt,:) = cell2mat(TAB(2:3,7));
        end
        
        for nperm=1:npermutation
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... sample %4.0f perm %4.0f',nt,nperm)
            pdm = data(all_perms(:,nperm),:);
            [~,TAB,~] = anovan(pdm(:,nt), group,'model','full','display','off','continuous',contVarIndex);
            if strcmp(model,'full')
                pd(nperm,nt,:) = cell2mat(TAB(2:4,6));
%                 pdpv(nperm,nt,:) = cell2mat(TAB(2:4,7));
            else
                pd(nperm,nt,:) = cell2mat(TAB(2:3,6));
%                 pdpv(nperm,nt,:) = cell2mat(TAB(2:3,7));
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


