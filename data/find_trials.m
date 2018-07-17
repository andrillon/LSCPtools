function idx_trials=find_trials(myconditions,mystring)

idx_trials=find(~(cellfun(@isempty,regexp(myconditions,mystring))));