%% Hyperspectral unmixing via L1_2 Sparsity-constrained 
% non-negative matrix factorization.
%-----------------------------------------------------------------------------------
    % Paper:
    % Hyperspectral Unmixing Via L1_2 Sparsity-constrained 
    % Nonnegative Matrix Factorization
%-----------------------------------------------------------------------------------
function [S, A, ARc] = hyperNmfASCL1_2(X, AInit, SInit, tolObj, maxIter, fDelta)
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
% Currenly, we omit this operation since we already know the number of 
% endmember in synthetic data.
emNum = 3;


% Estimate the weight parameter fLamda according to the sparsity measure
% over X.
bandNum = size(X, 1);
sampleNum = size(X, 2);
sqrtSampleNum = sqrt(sampleNum);
tmp = 0;
for l=1:bandNum
    xl =  X(l, :);
    tmp = tmp + ( sqrtSampleNum - (norm(xl,1)/norm(xl,2)) ) / ( sqrtSampleNum -1 );
end
fLamda = tmp / sqrt(bandNum);


% Initialize A and S by randomly selecting entries in the interval [0 1].
% Rescale each column of S to unit norm.
ARecord = zeros(emNum, maxIter, bandNum);
A = AInit;
S = SInit;
iterNum = 1;


% Run iterations.
Xf = [ X; fDelta*ones(1, sampleNum) ];
Af = [ A; fDelta * ones(1, emNum) ];
newObj = 0.5 * norm( (Xf - Af*S), 2 )^2 + fLamda * fNorm(S, 1/2);
oldObj = 0;
dispStr = ['Iteration ' num2str(iterNum),...
           ' loss = ' num2str(newObj)];
disp(dispStr);

for i = 1:emNum
    ARecord(i, iterNum, :) = A(1:bandNum, i);
end

while ( (oldObj-newObj)^2 >tolObj && (iterNum < maxIter) )
    oldObj = newObj;
    % update A
    A = Af .* (Xf*S') ./ (Af*S*S');
    S1_2 = S.^(-1/2);
    upperLimit = 10;
    S1_2(find(S1_2>upperLimit)) = upperLimit;
%     S = S .* (A'*Xf) ./ (A'*A*S);
    S = S .* (A'*Xf) ./ (A'*A*S + 0.5*fLamda*S1_2);
    Xf = [ X; fDelta*ones(1, sampleNum) ];
    Af = [ A(1:bandNum,:); fDelta * ones(1, emNum) ];
    newObj = 0.5 * norm( (Xf - Af*S), 2 )^2 + fLamda * fNorm(S, 1/2);
    
    iterNum = iterNum + 1;
    dispStr = ['Iteration ' num2str(iterNum),...
                ' loss = ' num2str(newObj)];
    disp(dispStr);
    
    for i = 1:emNum
        ARecord(i, iterNum, :) = A(1:bandNum, i);
    end
end

ARc = ARecord(:, 1:iterNum, :);

end