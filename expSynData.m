%% Exprimenet on synthetic data.

% Generate data.
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 5;
emNum = 6;
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
vca_ = false;
if vca_
    HVca = hyperVca(V', emNum);
end

% Run Nfindr
indxNfindr = nFindr(V, emNum);
HNfindr = V(indxNfindr, :);
sadNfindr = sadEms(HTrue, HNfindr, emNum);
sadSet = sadNfindr;
legendBar = { 'N-Findr' };

% Run nmfAbundance to find initial well-conditioned abundance.
HI = HNfindr;
alpha = 1;
tol = 0.1;
maxIter = 5000;
[WI, EI] = nmfAbundance(V, emNum, HI,...
                        alpha, tol, maxIter);

% Run hyperNmfMDC.
mdc_ = true;
if mdc_
    [ WMdc, HMdc, HRcMDC, EMDC] = ...
        hyperNmfMDC(...
        V, emNum, WI, HI, ...
        0.0002, 4,...
        0.01, 30000 );
    sadMdc = sadEms(HTrue, HMdc, emNum);
    sadSet = [sadSet, sadMdc];
    legendBar{ length(legendBar)+1 } = 'MDC';
end

% Run hyperNmfMVC
mvc_ = true;
if mvc_
    [UU, SS, WW] = svd(V');
    prinComp = pca(V);
    meanData = mean(HTrue, 1);

    [HMvc, WMvc] = hyperNmfMVC(V', HI', WI', ...
        HI', UU, prinComp, meanData, ...
        0.015, 0.001, 100, ...
        0, 2, 1);
    sadMvc = sadEms(HTrue, HMvc', emNum);
    sadSet = [sadSet, sadMvc];
    legendBar{ length(legendBar)+1 } = 'ICE';
end

% Run hyperNmfASCL1_2
spl1_2_ = false;
if spl1_2_
    [ WASCL1_2, HASCL1_2, HRcL1_2, errRcL1_2, objRcL1] = ...
        hyperNmfASCL1_2(...
        V', HI', WI',...
        0.001,... % tolObj
        20000,... % maxIter
        20 ... %fDelta
        );
    sadASCL1_2 = sadEms(HTrue, HASCL1_2, emNum);
    sadSet = [sadSet, sadASCL1_2];
    legendBar{ length(legendBar) } = 'ASCL_{1/2}';
end

% Run hyperNmfASCL1
spl1_ = false;
if spl1_
    [ WASCL1, HASCL1, HRcL1, errRcL1, objRcL1] = ...
        hyperNmfASCL1(...
        V', HI', WI',...
        0.001,... % tolObj
        20000,... % maxIter
        20 ... %fDelta
        );
    sadASCL1 = sadEms(HASCL1, HASCL1_2, emNum);
    sadSet = [sadSet, sadASCL1];
    legendBar{ length(legendBar) } = 'ASCL_1';
end

% Visualize result
figure;
hold on
% true data
scatter(V(:,1), V(:,2), 'm+');
hTrue = scatter(HTrue(:, 1), HTrue(:, 2), 100, 'filled', 'r');
handlesScatter = hTrue;
legendScatter = {'True'};

% VCA
if vca_
    hVCA = scatter(HVca(1,:), HVca(2,:), 100, 'filled', 'o');
    handlesScatter = [handlesScatter, hVCA];
    legendScatter{ length(legendScatter)+1 } = 'VCA';
end

% Nfindr
hNfindr = scatter(HNfindr(:, 1), HNfindr(:, 2), 100, 'filled', 'm');
handlesScatter = [handlesScatter, hNfindr];
legendScatter{ length(legendScatter)+1 } = 'N-Findr';

% MDC
if mdc_
    bandIndx1 = 1;
    bandIndx2 = 2;
    hMDC = scatter( HMdc(:,bandIndx1), HMdc(:,bandIndx2) , 100, 'filled','k');
    for i=1:bandNum
        plot(HRcMDC(i, :,bandIndx1), HRcMDC(i, :,bandIndx2), 'm-', 'MarkerSize', 5);
    end
    VMdc = WMdc * HMdc;
    scatter(VMdc(:,1), VMdc(:,2), 'o');
    handlesScatter = [handlesScatter, hMDC];
    legendScatter{ length(legendScatter)+1 } = 'MDC';
end

% MVC
if mvc_
    hMVC = scatter( HMvc(bandIndx1,:), HMvc(bandIndx2,:) , 100, 'filled','b');
    handlesScatter = [handlesScatter, hMVC];
    legendScatter{ length(legendScatter)+1 } = 'ICE';
end

% ASCL1_2
if spl1_2_
    hASCL1_2 = scatter( HASCL1_2(bandIndx1, :), HASCL1_2(bandIndx2,:) , 100, 'filled','c');
    for i = 1:emNum
        plot(HRcL1_2(i, :,bandIndx1), HRcL1_2(i, :,bandIndx2), 'r-.', 'MarkerSize', 5);
    end
    handlesScatter = [handlesScatter, hASCL1_2];
    legendScatter{ length(legendScatter)+1 } = 'ASCL_{1/2}';
end

% ASCL1
if spl1_
    hASCL1 = scatter( HASCL1(bandIndx1, :), HASCL1(bandIndx2,:) , 100, 'filled','g');
    for i = 1:emNum
        plot(HRcL1(i, :,bandIndx1), HRcL1(i, :,bandIndx2), 'r-.', 'MarkerSize', 5);
    end
    handlesScatter = [handlesScatter, hASCL1];
    legendScatter{ length(legendScatter)+1 } = 'ASCL_1';
end

legend(handlesScatter, legendScatter)


% Computer reconstruction error.
alpha = 1;
tol = 0.1;
maxIter = 5000;
figure;
hold on

hBar = bar(sadSet);
xlabel('endmember', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('SAD', 'FontSize', 15, 'FontWeight', 'bold')
legend(hBar, legendBar)
