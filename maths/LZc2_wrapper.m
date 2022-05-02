function [LZ, perm, HL_nrm] = LZc2_wrapper(X,savename,verbose)
% Wrapper for LZc2. Outputs Lempel-Ziv (LZ) complexity and the permuted
% complexity and entropy normalisation values.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - March 2016
%
% INPUTS:
%
% X        - EEG data (channels x samples x trials)
% savename - path and name to save temporary results to (optional).
% verbose  - Print progress? (default = 0)
%
% OUTPUTS:
%
% LZ       - Lempel-Ziv features.
% perm     - Features from randomly permuted EEG.
%       .c - LZ complexity.
%      .mw - LZ maximum no. of words of any word length.
%      .ml - LZ maximum word length.
% HL_nrm   - Entropy normalisation constatnt for LZc.

%% Preallocate and constants
if nargin < 3
    verbose = 0;
end

if nargin < 2 || (nargin>=2 && isempty(savename))
   dosave = 0; 
end

[D,N,Nwin] = size(X);
LZ.c = NaN(Nwin,1);
LZ.mw = NaN(Nwin,1);
LZ.ml = NaN(Nwin,1);
perm.c = NaN(Nwin,1);
perm.mw = NaN(Nwin,1);
perm.ml = NaN(Nwin,1);
HL_nrm = NaN(Nwin,1);
L = D*N;
const = log2(L)/L;

print_prc = 0.01;
print_ind = 0;


%% Loop over windows
tic1 = tic;
for n=1:Nwin
    Y = Pre(X(:,:,n));
    
    % binarising
    s = binarise(Y);
    
    % calculating LZc
    [LZ.c(n), LZ.mw(n), LZ.ml(n)] = LZc2(s);
    
    % calculating permuted LZc and entropy for normalisation.
    s_perm = s(randperm(length(s)));
    [perm.c(n), perm.mw(n), perm.ml(n)] = LZc2(s_perm);
    p1 = mean(s=='1');
    HL = -p1*log2(p1) - (1-p1)*log2(1-p1);
    HL_nrm(n) = const/HL;
    
    if n/Nwin > print_prc*print_ind
        print_ind = print_ind + 1;
        if verbose
            fprintf(['\n\n\n\n\n'...
                '%.1f %% done. Time spent: %.2f s. Estimated total time: %.0f s.\n\n\n\n\n']...
                ,n/Nwin*100,toc(tic1),toc(tic1)*Nwin/n)
        end
        if dosave
            save(savename,'LZ', 'perm', 'HL_nrm')
        end
    end
    
end
if dosave
    save(savename,'LZ', 'perm', 'HL_nrm')
end
end

function Y = Pre(X)
% Detrending data in X rowwise.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - October 2015
[D,N] = size(X);
Y = zeros(D,N);

for d = 1:D
    temp=detrend( X(d,:) - mean(X(d,:)) );
    Y(d,:) = temp;
end
end

function s = binarise(X)
% Binarises a multidimensional input to a single string of 1's and 0's
% using hilbert transform.
% Input: Continuous multidimensional time series.
% Output: One string being the binarized input matrix concatenated comlumn-
% by-column.
% Translated from python code supplied for [1].
% [1] Schartner et. al. Complexity of Multi-Dimensional Spontaneous EEG
% Decreases during Propofol Induced General Anaesthesia. Plos One (2015).
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - October 2015

[D,N] = size(X);
A = zeros(D,N);
for d = 1:D
    A(d,:) = abs(hilbert(X(d,:)));
end
thresh = mean(A,2);

s = char(zeros(D*N,1));
ind = 0;
for n = 1:N
    for d = 1:D
        ind = ind + 1;
        if A(d,n) > thresh(d)
            s(ind) = '1';
        else
            s(ind) = '0';
        end
    end
end

end