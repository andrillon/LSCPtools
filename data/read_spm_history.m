function read_spm_history(D)

hist=D.history;
for n=1:length(D.history)
    stepname=hist(n).fun;
    fprintf('step%g: %s\n',n,stepname)
    steparg=hist(n).args;
    if isfield(steparg,'filter')
        if length(steparg.filter.PHz)==1
        fprintf('\tband:%s\n\ttype:%s\n\torder:%g\n\tdir:%s\n\tHz:%g\n',...
            steparg.filter.band,steparg.filter.type,steparg.filter.order,steparg.filter.dir,steparg.filter.PHz(1))
        elseif length(steparg.filter.PHz)==2
             fprintf('\tband:%s\n\ttype:%s\n\torder:%g\n\tdir:%s\n\tHz:%g-%g\n',...
            steparg.filter.band,steparg.filter.type,steparg.filter.order,steparg.filter.dir,steparg.filter.PHz(1),steparg.filter.PHz(2))
        end
    elseif isfield(steparg,'montage')
        
    end
end