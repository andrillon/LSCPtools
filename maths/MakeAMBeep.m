function [beep,samplingRate] = MakeAMBeep(freqCarrier,freqAM,rateAM,duration,samplingRate)

ModRate=rateAM;

if nargin<3 || isempty(duration)
	error('Usage: beep=MakeAMBeep(freqCarrier,freqAM,duration,[samplingRate]);')
end
if nargin<4 || isempty(samplingRate)
	samplingRate = Snd('DefaultRate');
end

if length(duration)==1
    beep = sin(2*pi*freqCarrier*(0:duration*samplingRate)/samplingRate);
    
    amplMod = (1+ModRate*sin(2*pi*freqAM*(0:duration*samplingRate)/samplingRate));
    
    beep=beep.*amplMod;
elseif length(duration)==2
    totduration=sum(duration);
    
    beep = sin(2*pi*freqCarrier*(0:totduration*samplingRate)/samplingRate);
    
    amplMod = (1+ModRate*sin(2*pi*freqAM*(0:totduration*samplingRate)/samplingRate));
    
    beep=beep.*amplMod;
    
    ramp=[0:1/(duration(1)*samplingRate):1 ones(1,length(beep)-length(0:1/(duration(1)*samplingRate):1))];
    
    beep=beep.*ramp;
    
elseif length(duration)==3
    
    totduration=sum(duration);
    
    beep = sin(2*pi*freqCarrier*(0:totduration*samplingRate)/samplingRate);
    
    amplMod = (1+ModRate*sin(2*pi*freqAM*(0:totduration*samplingRate)/samplingRate));
    
    beep=beep.*amplMod;
    
    ramp=[0:1/(duration(1)*samplingRate):1 ones(1,length(beep)-length(0:1/(duration(1)*samplingRate):1)-length(0:1/(duration(3)*samplingRate):1)) 1:-1/(duration(3)*samplingRate):0];
    
    beep=beep.*ramp;
end


% sf = 22050;                        % sample frequency (Hz)
% d = 1.0;                           % durati0n (s)
% n = sf * d;                        % number of samples
% 
% % set carrier
% cf = 1000;                         % carrier frequency (Hz)
% c = (1:n) / sf;                    % carrier data preparation
% c = sin(2 * pi * cf * c);          % sinusoidal modulation
% 
% % set modulator
% mf = 5;                            % modulator frequency (Hz)
% mi = 0.5;                          % modulator index
% m = (1:n) / sf;                    % modulator data preparation
% m = 1 + mi * sin(2 * pi * mf * m); % sinusoidal modulation
% 
% % amplitude modulation
% s = m .* c;                        % amplitude modulation
% 
% % sound presentation
% sound(s, sf);                      % sound presentation
% pause(d + 0.5);                    % waiting for sound end
% plot(s);