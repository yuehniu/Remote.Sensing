%% Exprimenet on synthetic data.

% Generate data.
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 4;
emNum = 4;
HTrue = zeros(emNum, bandNum);
HI = zeros(emNum, bandNum);
HMdc = zeros(emNum, bandNum);
WTrue = zeros(dataSize, emNum);
WI = zeros(dataSize, emNum);
WMdc = zeros(dataSize, emNum);
HTrue = abs( randn( emNum, bandNum ) );
[V, W_true] = create4(dataSize, HTrue);

figure;
scatter(V(:,1), V(:,2));
xlabel('band 1');
ylabel('band 2');

% Run VCA
HVca = hyperVca(V', emNum);

% Run Nfindr
indxNfindr = nFindr(V, emNum);
HNfindr = V(indxNfindr, :);

% Run nmfAbundance to find initial well-conditioned abundance.
HI = HNfindr;
alpha = 1;
tol = 0.1;
maxIter = 5000;
[WI, EI] = nmfAbundance(V, emNum, HI,...
                        alpha, tol, maxIter);

% Run hyperNmfMDC.
[ WMdc, HMdc, HRcMDC, EMDC] = ...
    hyperNmfMDC(...
    V, emNum, WI, HI, ...
    0.001, 4,...
    0.01, 30000 );

% Run hyperNmfMVC
[UU, SS, WW] = svd(V');
prinComp = pca(V);
meanData = mean(HTrue, 1);

[HMvc, WMvc] = hyperNmfMVC(V', HI', WI', ...
    HI', UU, prinComp, meanData, ...
    0.015, 0.001, 30000, ...
    0, 2, 1);

% Run hyperNmfASCL1_2
[ WASCL1_2, HASCL1_2, HRcL1_2, errRcL1_2, objRcL1] = ...
    hyperNmfASCL1_2(...
    V', HI', WI',...
    0.001,... % tolObj
    20000,... % maxIter
    20 ... %fDelta
    );

% Run hyperNmfASCL1
[ WASCL1, HASCL1, HRcL1, errRcL1, objRcL1] = ...
    hyperNmfASCL1(...
    V', HI', WI',...
    0.001,... % tolObj
    20000,... % maxIter
    20 ... %fDelta
    );


% Visualize result
figure;
hold on
% true data
scatter(V(:,1), V(:,2));
hTrue = scatter(HTrue(:, 1), HTrue(:, 2), 'filled', 'r');
% VCA
hVCA = scatter(HVca(1,:), HVca(2,:), 100, 'filled', 'o');
% Nfindr
hNfindr = scatter(HNfindr(:, 1), HNfindr(:, 2), 'filled', 'm');
% MDC
bandIndx1 = 1;
bandIndx2 = 2;
hMDC = scatter( HMdc(:,bandIndx1), HMdc(:,bandIndx2) , 'filled','k');
for i=1:bandNum
    plot(HRcMDC(i, :,bandIndx1), HRcMDC(i, :,bandIndx2), 'm-', 'MarkerSize', 5);
end
% MVC
hMVC = scatter( HMvc(bandIndx1,:), HMvc(bandIndx2,:) , 'filled','b');
% ASCL1_2
hASCL1_2 = scatter( HASCL1_2(bandIndx1, :), HASCL1_2(bandIndx2,:) , 'filled','c');
for i = 1:emNum
    plot(HRcL1_2(i, :,bandIndx1), HRcL1_2(i, :,bandIndx2), 'r-.', 'MarkerSize', 5);
end
% ASCL1
hASCL1 = scatter( HASCL1(bandIndx1, :), HASCL1(bandIndx2,:) , 'filled','g');
for i = 1:emNum
    plot(HRcL1(i, :,bandIndx1), HRcL1(i, :,bandIndx2), 'r-.', 'MarkerSize', 5);
end

legend([hTrue,hVCA,hNfindr,hMDC,hMVC,hASCL1_2,hASCL1], 'True', 'VCA', 'Nfindr', 'MDC', 'MVC', 'ASCL1/2', 'ASCL1')


% Computer reconstruction error.
alpha = 1;
tol = 0.1;
maxIter = 5000;
WVca = nmfAbundance(V, emNum, HVca',...
                        alpha, tol, maxIter);
errVca = sum(sum((V-WVca*HVca').^2));
errNfindr = sum(sum((V-WI*HNfindr).^2));
errMdc = sum(sum((V-WMdc*HMdc).^2));
errMvc = sum(sum((V-WMvc'*HMvc').^2));
errASCL1_2 = sum(sum((V-WASCL1_2'*HASCL1_2(1:emNum,:)').^2));
errASCL1 = sum(sum((V-WASCL1'*HASCL1(1:emNum,:)').^2));
figure;
hold on
heVCA = plot(1, errVca, '*', 'MarkerSize', 10, 'Color', 'r');
heNfindr = plot(1, errNfindr, '*', 'MarkerSize', 10, 'Color', 'm');
heMdc = plot(1, errMdc, '*', 'MarkerSize', 10, 'Color', 'k');
heMvc = plot(1, errMvc, '*', 'MarkerSize', 10, 'Color', 'b');
heASCL1_2 = plot(1, errASCL1_2, '*', 'MarkerSize', 10, 'Color', 'c');
heASCL1 = plot(1, errASCL1, '*', 'MarkerSize', 10, 'Color', 'g'); 
legend([heVCA,heNfindr,heMdc,heMvc,heASCL1_2,heASCL1], 'VCA', 'Nfindr', 'MDC', 'MVC', 'ASCL1/2', 'ASCL1');
