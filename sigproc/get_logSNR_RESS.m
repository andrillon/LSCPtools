function [snrR, snrE, hz, map2plot, ress_ts1]=get_logSNR_RESS(data,SR,param)
% Function that extracts the logSNR of the data (frequency tagging
% experiments) using the RESS method (Cohen % Gulbinaite Neuroimage 2016)
%
% data: EEG/ECoG data (channel x time x trial)
% SR: sampling rate in Hz
% param: structure with fields:
%   - peakfreq1: Target Frequency
%   - peakwidt: FWHM at peak frequency
%   - neighfreq: distance of neighboring frequencies away from peak frequency, +/- in Hz
%   - neighwidt: FWHM of the neighboring frequencies
%   - fft_res: Resolution of FFT (e.g. 0.1 Hz)
%   - mindist: Min distance between peaks (0.5Hz)
%   - lowerfreq: Lower frequency bound (2Hz)

peakfreq1=param.peakfreq1;
peakwidt=param.peakwidt;
neighfreq=param.neighfreq;
neighwidt=param.neighwidt;
fft_res=param.fft_res; %0.1 default
mindist=param.mindist; %0.5Hz default
lowerfreq=param.lowerfreq; %2Hz default


%%% Start RESS for a given condition
% FFT parameters
nfft = ceil(SR/fft_res ); % .1 Hz resolution

% extract EEG data
dataX = mean(abs(fft(data(:,:,:),nfft,2)/size(data,2)).^2,3);
hz    = linspace(0,SR,nfft);

% compute covariance matrix at peak frequency
fdatAt = filterFGx(data,SR,peakfreq1,peakwidt);
fdatAt = reshape( fdatAt(:,:,:), size(data,1),[] );
fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
covAt  = (fdatAt*fdatAt')/size(data,2);

% compute covariance matrix for lower neighbor
fdatLo = filterFGx(data,SR,peakfreq1+neighfreq,neighwidt);
fdatLo = reshape( fdatLo(:,:,:), size(data,1),[] );
fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
covLo  = (fdatLo*fdatLo')/size(data,2);

% compute covariance matrix for upper neighbor
fdatHi = filterFGx(data,SR,peakfreq1-neighfreq,neighwidt);
fdatHi = reshape( fdatHi(:,:,:), size(data,1),[] );
fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
covHi  = (fdatHi*fdatHi')/size(data,2);

% perform generalized eigendecomposition. This is the meat & potatos of RESS
[evecs,evals] = eig(covAt,(covHi+covLo)/2);
[~,comp2plot] = max(diag(evals)); % find maximum component
evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)

% extract components and force sign
% maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
maps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
[~,idx] = max(abs(maps(:,comp2plot))); % find biggest component
maps = maps * sign(maps(idx,comp2plot)); % force to positive sign
map2plot = maps(:,comp2plot);

% reconstruct RESS component time series
ress_ts1 = zeros(size(data,2),size(data,3));
for ti=1:size(data,3)
    ress_ts1(:,ti) = evecs(:,comp2plot)'*squeeze(data(:,:,ti));
end

%%% compute SNR spectrum for all channels
ressx = mean(abs( fft(ress_ts1(:,:),nfft,1)/size(data,2)).^2,2);


[snrR,snrE] = deal(zeros(size(hz)));
skipbins =  mindist/fft_res; % .5 Hz, hard-coded!
numbins  = lowerfreq/fft_res+skipbins; %  2 Hz, also hard-coded!

% loop over frequencies and compute SNR
for hzi=numbins+1:length(hz)-numbins-1
    numer = ressx(hzi);
    denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
    snrR(1,hzi) = numer./denom;
    
    for nch=1:size(data,1)
        elecx = dataX(nch,:,:);
        numer = elecx(hzi);
        denom = mean( elecx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snrE(nch,hzi) = numer./denom;
    end
end
