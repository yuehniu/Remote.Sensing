%% Hyperspectral unmixing via L1_2 Sparsity-constrained 
% non-negative matrix factorization.
%-----------------------------------------------------------------------------------
    % Paper:
    % Hyperspectral Unmixing Via L1_2 Sparsity-constrained 
    % Nonnegative Matrix Factorization
%-----------------------------------------------------------------------------------
function [S, A, ARc, errRc, objRc] = hyperNmfASCL1_2(X, AInit, SInit, tolObj, maxIter, fDelta)
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
%     ARc: iteraive record for endmember matrix (emNum * iterNum * bandNum)
%     errRc: iterative record for error (iterNum).
%     objRc: iterative record for object value (iterNum).

% Estimate the number of enNum endmembers using the HySime algorithm.
% Currenly, we omit this operation since we already know the number of 
% endmember in synthetic data.
emNum = size(AInit, 2);


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

% Record iteration.
errRc = zeros(1, maxIter);
objRc = zeros(1, maxIter);
ARecord = zeros(emNum, maxIter, bandNum);

% fLamda should be rescale to the level of spectral sample value
fLamda = fLamda / 500;


% Initialize A and S by randomly selecting entries in the interval [0 1].
% Rescale each column of S to unit norm.
A = AInit;
S = SInit;
iterNum = 1;


% Run iterations.
Xf = [ X; fDelta*ones(1, sampleNum) ];
Af = [ A; fDelta * ones(1, emNum) ];
err = 0.5 * norm( (Xf(1:bandNum,:) - Af(1:bandNum,:)*S), 2 )^2;
newObj = err + fLamda * fNorm(S, 1);
oldObj = 0;
dispStr = ['Iteration ' num2str(iterNum),...
           ' loss = ' num2str(newObj)];
disp(dispStr);

% record iteration.
errRc(iterNum) = err;
objRc(iterNum) = newObj;
for i = 1:emNum
    ARecord(i, iterNum, :) = A(1:bandNum, i);
end

while ( err >tolObj && (iterNum < maxIter) )
    oldObj = newObj;
    % update A
    A = Af .* (Xf*S') ./ (Af*S*S');
    
    % update S
    % S = S .* (A'*Xf) ./ (A'*A*S);
    S = S .* (A'*Xf) ./ (A'*A*S + fLamda);
    % Xf = [ X; fDelta*ones(1, sampleNum) ];
    % Af = [ A(1:bandNum,:); fDelta * ones(1, emNum) ];
    Af = A;
    err = 0.5 * norm( (Xf(1:bandNum,:) - Af(1:bandNum,:)*S), 2 )^2;
    newObj = err + fLamda * fNorm(S, 1);
    
    iterNum = iterNum + 1;
    dispStr = ['Iteration ' num2str(iterNum),...
                ' loss = ' num2str(newObj)];
    disp(dispStr);
    
    % record iteration.
    errRc(iterNum) = err;
    objRc(iterNum) = newObj;
    for i = 1:emNum
        ARecord(i, iterNum, :) = A(1:bandNum, i);
    end
end

ARc = ARecord(:, 1:iterNum, :);
errRc = errRc(1, 1:iterNum);
objRc = objRc(1, 1:iterNum);

end