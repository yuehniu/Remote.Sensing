%% Visualize simplex and discover some properties.

% Generate a simplex
endNum = 3;
sampleNum = 6400;
A = rand(endNum, sampleNum);
A = A ./ repmat(sum(A, 1), endNum, 1);

% Scatter A
showFlag = true;
if showFlag
    figure;
    indx1 = 1;
    indx2 = 2;
    scatter(A(indx1, :), A(indx2, :))
end


% Generate spectral samples
bandNum = 2;
M = rand(bandNum, endNum);
X = M * A;

% scatter X
showFlag = true;
if showFlag
    figure;
    indx1 = 1;
    indx2 = 2;
    scatter(X(indx1, :), X(indx2, :))
end


% Generate simplex cone
gamma = 1:1:100;
figure;
hold on
for i=1:length(gamma)
    Xgamma = gamma(i) * X;
    indx1 = 1;
    indx2 = 2;
    scatter(Xgamma(indx1, :), Xgamma(indx2, :))
end