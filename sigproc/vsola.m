function op = vsola(ip, tsm_factor, P)
% Implementation of the variable parameter  synchronised overlap add VSOLA algorithm 
% This implementation makes use of the standard SOLA algorithm (Rocus and Wilgus, ICASSP 1986) and some efficient paramter settings for the SOLA algorithm (Dorran, Lawlor and Coyle, ICASSP 2003) and (Dorran, Lawlor and Coyle, DAFX 2003)
% Given an input, ip, the longest likely pitch period, P, and  and a time scale modification factor, tsm_factor, return an output, op, that is a time scaled version of the input. The synthesis_overlap_size is the length, in samples of the lowest likely fundametal period of the input.
% for speech synthesis_overlap_size is set to (16/1000)*fs samples and for music synthesis_overlap_size is typically set to (20/1000)*fs samples 
%
% David Dorran, Audio Research Group, Dublin Institute of Technology
% david.dorran@dit.ie
% http://eleceng.dit.ie/dorran
% http://eleceng.dit.ie/arg
%
 
% make sure input is mono and transpose if necessary
[r, c] = size(ip);
if r > 1
    ip = ip';
end;    
[r, c] = size(ip);
if r > 1
    disp('Note :only works on mono signals');
    op = [];
    return
end;
 
% initialize the values of analysis_frame_length, analysis_window_offset, synthesis_frame_offset and length_of_overlap
desired_tsm_len = round(length(ip)*tsm_factor);
P = round(P); %longest likely pitch period in samples
Lmax = round(P * 1.5);% found this to a reasonable value for the Lmax- Lmax is the duration over which the correlation function is applied
stationary_length = (P * 1.5); % found this to a reasonable value for the stationary length - This is the max duration that could be discarded/repeated
 
analysis_window_offset = round((stationary_length - P)/abs(tsm_factor - 1)); % this equation was derived.
synthesis_window_offset = round(analysis_window_offset * tsm_factor);
analysis_frame_length = round(Lmax + synthesis_window_offset);
number_of_analysis_frames = floor((length(ip)- analysis_frame_length)/analysis_window_offset);
 
if number_of_analysis_frames < 2 %not much time-scaling being done just return the input
    op = ip;
    return;
end;
 
%the next two lines just ensure that the last frame finishes at the very end of the signal (not essential)
zpad = zeros(1, (number_of_analysis_frames*analysis_window_offset) + analysis_frame_length - length(ip));
ip = [ip zpad];
 
%initialize the output
op = zeros(1, desired_tsm_len);
%initialize the first output frame
op(1 : analysis_frame_length) = ip(1 : analysis_frame_length);
 
min_overlap = round(Lmax - P); %ensure that there is some minimum overlap
count = 0;
 
% Loop for the 2nd analysis frame to the number_of_analysis_frames
for m = 1 : number_of_analysis_frames
     
    %grab the mth input frame
    ip_frame = ip(analysis_window_offset * m : (analysis_window_offset * m) + analysis_frame_length - 1);
     
    %grab the maximum overlapping segments from the inout frame and the current output
    seg_1 = op(round(synthesis_window_offset*(m-1))+analysis_frame_length - Lmax : round(synthesis_window_offset*(m-1))+analysis_frame_length -1);
    seg_2 = ip_frame(1: Lmax);
     
    %compute the correlation of these segments
    correlation   = xcorr(seg_2, seg_1,'coeff');
 
    %Find the best point to overlap (opt_overlap_length) making sure not to exceed the maximum or go below the minimum overlap.
    correlation(length(correlation) - Lmax -1: length(correlation)) = -100;
    correlation(1: min_overlap) = -100;    
    [max_correlation, opt_overlap_length] = max(correlation);
     
    if(max_correlation == 0)
        opt_overlap_length = Lmax;
    end;
%     offset = Lmax - opt_overlap_length;
%     if ((offset + analysis_window_offset -  synthesis_window_offset) >= 0 & (offset + analysis_window_offset -  synthesis_window_offset) <= P)
%         count = count +1;
%     end;
     
    % append mth analysis frame to the current synthesised output using a linear cross fade
    ov_seg = linear_cross_fade(seg_1, ip_frame, opt_overlap_length);
    ov_len =(round(synthesis_window_offset*m)+analysis_frame_length) - (round(synthesis_window_offset*(m-1))+analysis_frame_length - Lmax) + 1;
    ov_seg = ov_seg(1:ov_len);
    op(round(synthesis_window_offset*(m-1))+analysis_frame_length - Lmax: round(synthesis_window_offset*m)+analysis_frame_length) = ov_seg;
     
end; % end of for loop
 
% linear cross fade the first segment with the second segment given a certain amount of overlap
% |----------seg1----------|
%              |---------------seg2-----------------------|
%              |--overlap-|
 
function op = linear_cross_fade(seg_1, seg_2, overlap_length, cross_fade_duration)
 
    error(nargchk(3,4,nargin));
    if nargin < 4
        cross_fade_duration = overlap_length;
    end
    if cross_fade_duration > overlap_length
        cross_fade_duration = overlap_length;
    end
    % overlap the end of seg_1 with the start of seg_2 using a linear cross-fade
    if (length(seg_1) < overlap_length),
        seg_2 = seg_2(overlap_length - length(seg_1) + 1: length(seg_2));
        overlap_length = length(seg_1);
    end; % end of if statement
 
    if (length(seg_2) < overlap_length),
        seg_1 = seg_1(length(seg_1) - (overlap_length - overlap_length): length(seg_1));
        overlap_length = length(seg_2);
    end; % end of if statement
 
     
    overlapping_region = zeros(1, cross_fade_duration); % initialize the overlapping region
    seg_1 = seg_1(1: length(seg_1) - (overlap_length - cross_fade_duration));
     
    op = zeros(1, length(seg_1) + length(seg_2) - cross_fade_duration);
 
    if(overlap_length ~= 1)
        linear_ramp1 = (cross_fade_duration-1:-1:0) ./(cross_fade_duration-1);
        linear_ramp2 = (0:1:cross_fade_duration-1) ./(cross_fade_duration-1);
        overlapping_region = (seg_1(length(seg_1)-cross_fade_duration+1:length(seg_1)).*linear_ramp1) + (seg_2(1: cross_fade_duration).*linear_ramp2);
    end;
    op = [seg_1(1:length(seg_1)-cross_fade_duration) ,overlapping_region , seg_2(cross_fade_duration+1:length(seg_2))];
% END of linear_cross_fade function