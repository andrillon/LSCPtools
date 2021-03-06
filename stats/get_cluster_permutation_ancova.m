function [realpos realneg]=get_cluster_permutation_ancova(data,group,covar,contidx,montecarloalpha,clusteralpha,npermutation,sTime)

% Input
% - data: cell array (1 cell per condition:  channels * times * subjects
% - montecarloalpha
% - clusteralpha
% - npermutation
% - sTime
% - averageelecFlag

%% Bootstrap
% fprintf('Calculating permutations...');
data=data(~isnan(group),:);
group=group(~isnan(group));

ntime=size(data,2);
nsuj=size(data,1);
df1=length(unique(group))-1;
df2=nsuj-1;
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
        [~, anovatab, ~] = anovan(rdm(:,nt),[group' covar],'model','linear','display','off','continuous',contidx);
%     [~,TAB,~] = anova1(rdm(:,nt), group,'off');
    rd(nt) = anovatab{2,6};
    
    for nperm=1:npermutation
        pdm = data(randperm(size(data,1)),:);
        [~, anovatab, ~] = anovan(pdm(:,nt),[group' covar],'model','linear','display','off','continuous',contidx);
%         [~,TAB,~] = anova1(pdm(:,nt), [group covar],'off');
        pd(nperm,nt) = anovatab{2,6};
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


fprintf('\n');

pmonte = realneg.pmonte;
goodc = find(pmonte < montecarloalpha);
contrast = linspace(.5, 1, length(goodc));
for i = 1:length(goodc)
    ic = goodc(i);
    samples = realneg.clusters == ic;
    cint = [min(sTime(samples)) max(sTime(samples))];
    [peakv,peaki] = min(rd(samples));
    cintsamples = find(samples);
    peakt = sTime(cintsamples(peaki));
    fprintf('\t neg | p-value : %0.4f | time :  %1.3f %1.3f [peak : %1.3f]; ... \n', pmonte(ic), cint, peakt);
    
end

fprintf('\n')

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


