function [lzc, maxwords, maxlength] = LZc2(s)
% LZC computes the complexity measure, LZc, of a given string, s,
% consisting of of 1's and 0's. Complexity is given as the length of the
% estimated codebook. Translated from python code supplied for [1]. The
% function has been optimised wrt. the lookups in the codebook, but
% otherwise the computations are done as in [1]. Two new outputs have been
% defined as the maximum number of words of the same length, and the
% length of the longest word in the codebook.
% [1] Schartner et. al. Complexity of Multi-Dimensional Spontaneous EEG
% Decreases during Propofol Induced General Anaesthesia. Plos One (2015).
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - March 2016

%% Initialising
Lw = ceil(sqrt(2*length(s))); % max word length possible
d=cell(Lw,1); % codebook
d(:) = {''};
w =''; % the first word initialised as empty
lzc=0;

%% Iterating through each sample in s
for i=1:length(s)
    % add the next sample in s to the current word
    c = s(i); 
    wc = [w c];
    word_length = length(wc);
        
    % checks if word exists in codebook
    if any(ismember(d{word_length},wc));%~isempty(d{word_length}) && ismember(wc,d{word_length},'rows');%%%;%~isempty(d{word_length}) &&   %%
        w = wc;
    else
        % if not, the word is added to the codebook, and a new word is
        % started using this sample and the next
        lzc = lzc + 1;
        d{word_length}{end+1} = wc;
        w = c;
    end
end
%% The complexity is given a number of entries in codebook
maxwords = max(cellfun('length',d));
maxlength = find(~cellfun('isempty',d),1,'last');
end