%% Hyperspectral unmixing via L1_2 Sparsity-constrained 
% non-negative matrix factorization.
%-----------------------------------------------------------------------------------
    % Paper:
    % Hyperspectral Unmixing Via L1_2 Sparsity-constrained 
    % Nonnegative Matrix Factorization
%-----------------------------------------------------------------------------------
function [S, A, ARc, errRc, objRc] = hyperNmfMdcAscl1_2(X, AInit, SInit, tolObj, maxIter, dDelta, fDelta)
% Input:
%     X: Hyperspectral data maxtrix (sampleSize * bandNum).
%     AInit: Endmember intial matrix (emNum * bandNum).
%     SInit: Abundance initial matrix (sampleSize * emNum).
%     tolObj: Stop condition for difference of object function between two iteration.
%     maxIter: Maximun number of iterations.
%     dDelta: Factor that control the strength of distance constraint.
%     fDelta: Factor that control the strength of sum-to-one constraint.
%
% Output:
%     A: resultant endmember matrix (emNum * bandNum).
%     S: resultant abundance matrix (sampleSize * emNum).
%     ARc: iteraive record for endmember matrix (iterNum * emNum * bandNum)
%     errRc: iterative record for error (iterNum).
%     objRc: iterative record for object value (iterNum).

% Estimate the number of enNum endmembers using the HySime algorithm.
% Currenly, we omit this operation since we already know the number of 
% endmember in synthetic data.
emNum = size(AInit, 1);


% Estimate the weight parameter fLamda according to the sparsity measure
% over X.
bandNum = size(X, 2);
sampleNum = size(X, 1);
sqrtSampleNum = sqrt(sampleNum);
tmp = 0;
for l=1:bandNum
    xl =  X(:, l);
    tmp = tmp + ( sqrtSampleNum - (norm(xl,1)/norm(xl,2)) ) / ( sqrtSampleNum -1 );
end
fLamda = tmp / sqrt(bandNum);

% Transform matrix in MDC
M = eye(emNum);
M(emNum, :) = [1, zeros(1, emNum-1)];
P = M' + M;
I = eye(emNum);
% Record iteration.
errRc = zeros(1, maxIter);
objRc = zeros(1, maxIter);
ARecord = zeros(maxIter, emNum, bandNum);

% fLamda should be rescale to the level of spectral sample value
fLamda = fLamda / 500;


% Initialize A and S by randomly selecting entries in the interval [0 1].
% Rescale each column of S to unit norm.
A = AInit;
S = SInit;
iterNum = 1;


% Run iterations.
Xf = [ X fDelta*ones(sampleNum, 1) ];
Af = [ A fDelta * ones(emNum, 1) ];
err = 0.5 * norm( (Xf(:, 1:bandNum) - S*Af(:, 1:bandNum)), 2 )^2;
newObj = err + fLamda * fNorm(S, 1/2);
oldObj = 0;
dispStr = ['Iteration ' num2str(iterNum),...
           ' loss = ' num2str(newObj)];
disp(dispStr);

% record iteration.
errRc(iterNum) = err;
objRc(iterNum) = newObj;
for i = 1:emNum
    ARecord(iterNum, i, :) = A(i, 1:bandNum);
end

while ( err >tolObj && (iterNum < maxIter) )
    oldObj = newObj;
    % update A
    A = Af .* ((S'*Xf)+(dDelta * P * Af)) ./ ((S'*S+(2*dDelta * I))*Af);
    
    % update S
    lowLimit = 0.01;
    S(find(S<lowLimit))  = lowLimit;
    S1_2 = S.^(-1/2);
    % upperLimit = 10;
    % S1_2(find(S1_2>upperLimit)) = upperLimit;
    % S = S .* (A'*Xf) ./ (A'*A*S);
    S = S .* (Xf*A') ./ (S*A*A' + 0.5*fLamda*S1_2);
    % Xf = [ X; fDelta*ones(1, sampleNum) ];
    % Af = [ A(1:bandNum,:); fDelta * ones(1, emNum) ];
    Af = A;
    err = 0.5 * norm( (Xf(:, 1:bandNum) - S*Af(:,1:bandNum)), 2 )^2;
    newObj = err + fLamda * fNorm(S, 1/2);
    
    iterNum = iterNum + 1;
    dispStr = ['Iteration ' num2str(iterNum),...
                ' loss = ' num2str(newObj)];
    disp(dispStr);
    
    % record iteration.
    errRc(iterNum) = err;
    objRc(iterNum) = newObj;
    for i = 1:emNum
        ARecord(iterNum, i, :) = A(i, 1:bandNum);
    end
end

ARc = ARecord(1:iterNum, :, :);
errRc = errRc(1, 1:iterNum);
objRc = objRc(1, 1:iterNum);

end