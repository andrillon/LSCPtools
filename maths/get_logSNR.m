function [logSNR, faxis, logpow]=get_logSNR(data,SR,param)
% Function that extracts the logSNR of the data (frequency tagging
% experiments)
%
% data: EEG/ECoG data (channel x time x trial)
% SR: sampling rate in Hz
% param: structure with fields:
%   - method: 'fft' or 'taper'
%   - mindist: minimal expected distance between peaks (in Hz)
%   - numTaper: number of tapers if method if taper
if nargin<3
    param.method='fft';
    param.mindist=1;
end

if strcmp(param.method,'fft')
    T=size(data,2)/SR;
    df = 1/T;
    fNQ = SR/2;
    numfreq = length((0:df:fNQ));
    
    % pre-allocation for FFT
    logpow=nan(size(data,1),numfreq,size(data,3));
    logSNR=nan(size(data,1),numfreq,size(data,3));
end
for nCh=1:size(data,1)
    for nTr=1:size(data,3)
        
        signal=squeeze(data(nCh,:,nTr));
        if strcmp(param.method,'fft')
            %%%% Apply FFT
            % get power
            pow = (abs(fft(signal)).^2)/length(signal);
            % first half of data without negative frequencies
            pow = pow(1:min(floor(length(signal)/2)+1,length(signal)));
            % define df and fNQ
            fNQ = SR/2;
            faxis = (0:df:fNQ);
        elseif strcmp(param.method,'taper')
            % Apply taper methods (Chronux toolbox)
            T=size(data,2)/SR;
            df = 1/T;
            subparam=[];
            subparam.tapers=[df T 2*T*df-param.numTaper];
            subparam.Fs=SR;
            [pow,faxis]=mtspectrumc(signal,subparam);
        end
        
        %%%% Take the log
        thispow=log(pow);
        
        %%%% Compute the convolution kernel
        %         bw=1/(size(data,2)/SR); % bandwidth=1/T where T is the duration
        %         of the trial bw is actually df
        mindist=param.mindist;
        length_kernel=length(-mindist+df:df:mindist-df);
        kernel=ones(1,length_kernel);
        kernel((-mindist+df:df:mindist-df)>=-df & (-mindist+df:df:mindist-df)<=df)=0;
        kernel=kernel/sum(kernel);
        
        %%% apply kernel
        convlogpow= conv(thispow, kernel, 'same');
        logSNR(nCh,:,nTr)=thispow-convlogpow; % equivalent to logSNR=log(pow(fo))-1/n*SUM(log(pow(fi)))
        logpow(nCh,:,nTr)=thispow;
    end
end