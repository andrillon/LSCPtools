function [sym, count, smi, wsmi] = smi_and_wsmi(data, cfg)
    
    chan_sel  = cfg.chan_sel;
    data_sel  = cfg.data_sel;
    taus   = cfg.taus;
    kernel = cfg.kernel;

    over_trials = cfg.over_trials;

%%% autocomplete the weights matrix if not provided    
    if ~isfield(cfg, 'weights')
        nsymbols = factorial(kernel);
        wts = 1- (eye(nsymbols) | fliplr(eye(nsymbols)));
    else
        wts = cfg.weights;
    end

    %%% filtering and processing
    for tau=1:size(taus,2)
        tic
        filter_freq = cfg.sf/cfg.kernel/taus(tau);
        disp(['Filtering @ ' num2str(filter_freq) ' Hz'])
        ntrials = size(data,3);
        fdata = zeros(size(data));

        for trial=1:ntrials
            fdata(:,:,trial) = ft_preproc_lowpassfilter(squeeze(data(:,:,trial)),cfg.sf,filter_freq );  
        end
        disp(['Filtering done'])
        [sym{tau}, count{tau}, smi{tau}, wsmi{tau}] = smi_and_wsmi_mex(fdata(chan_sel, data_sel, :), taus(tau), kernel, wts, over_trials);
        toc
    end


end

