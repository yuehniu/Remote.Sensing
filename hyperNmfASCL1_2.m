%% Hyperspectral unmixing via L1_2 Sparsity-constrained 
% non-negative matrix factorization.
%-----------------------------------------------------------------------------------
    % Paper:
    % Hyperspectral Unmixing Via L1_2 Sparsity-constrained 
    % Nonnegative Matrix Factorization
%-----------------------------------------------------------------------------------
function [A, S] = hyperNmfASCL1_2(X, AInit, SInit, tolObj, maxIter, fDelta)
% Input:
%     X: Hyperspectral data maxtrix (bandNum * sampleSize).
%     AInit: Endmember intial matrix (bandNum * emNum).
%     SInit: Abundance initial matrix (emNum * sampleSize).
%     tolObj: Stop condition for difference of object function between two iteration.
%     maxIter: Maximun number of iterations.
%     fDelta: Factor that control the strength of sum-to-one constraint.
%
% Output:
%     A: resultant endmember matrix (bandNum * emNum).
%     S: resultant abundance matrix (emNum * sampleSize).


% Estimate the number of enNum endmembers using the HySime algorithm.


% Estimate the weight parameter fLamda according to the sparsity measure
% over X.


% Initialize A and S by randomly selecting entries in the interval [0 1].
% Rescale each column of S to unit norm.


% Run iterations.



end