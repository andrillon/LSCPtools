function BP = bandpass(timecourse, SamplingRate, low_cut, high_cut, filterOrder)
% BP = bandpass(timecourse, SamplingRate, low_cut, high_cut, filterOrder)

if (nargin < 5)
    filterOrder = 2;
end

[b, a] = butter(filterOrder, [(low_cut/SamplingRate)*2 (high_cut/SamplingRate)*2]);
BP = filtfilt(b, a, timecourse );

 
