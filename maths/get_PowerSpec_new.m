function [faxis, pow] = get_PowerSpec_new(signal, SamplingRate, lengthSec, DecibelsFlag ,plotFlag)
% For example for a 10sec segment sampled at millisec resolution use:
% [faxis, pow] = get_PowerSpec_new(mySignal, 1/1000, 10);
if (nargin < 5), plotFlag = 0; end
if (nargin < 4), DecibelsFlag = 0; end
% get power
pow = (abs(fft(signal)).^2)/length(signal);
% convert to decibels
if DecibelsFlag==1
    pow = 10*log10(pow/max(pow));
end
% first half of data without negative frequencies
pow = pow(1:min(floor(length(signal)/2)+1,length(signal)));
% define df and fNQ
df = 1/lengthSec;
fNQ = 1/SamplingRate/2;

faxis = (0:df:fNQ);
if (plotFlag)
    plot(faxis, pow);
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end