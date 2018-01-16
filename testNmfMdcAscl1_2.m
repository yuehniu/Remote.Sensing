%% nmf_MDC_test with more than 3 endmembers
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 4;
endNum = 4;
HTrue = zeros(endNum, bandNum);
HI = zeros(endNum, bandNum);
HMdcAscl1_2 = zeros(endNum, bandNum);
WTrue = zeros(dataSize, endNum);
WI = zeros(dataSize, endNum);
WMdcAscl1_2 = zeros(dataSize, endNum);
% for i = 1:exp_times
    % generate true W, H, V
    HTrue = abs( randn( endNum, bandNum ) );
    [V, W_true] = create4(dataSize, HTrue);
%     max_V = max(V);
%     for i = 1:bandNum
%        V(:,i) = V(:,i) / (max_V(i));
%     end
    figure;
    scatter(V(:,1), V(:,2));
    xlabel('band 1');
    ylabel('band 2');
    xlim([0,1])
    ylim([0,1])
    
    %% find initial H using n_findr
    HInitIndx = nFindr(V, endNum);
    HI = V(HInitIndx, :);
    
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;
    [WI, E_I] = nmfAbundance(V, endNum, HI,...
                        alpha, tol, maxIter);
    VNmf = WI * HI;
    
    figure;
    scatter(VNmf(:,1), VNmf(:,2), 'c', 'full'); hold on
    scatter(V(:,1), V(:,2), 'r')
    scatter(HI(:, 1), HI(:, 2), 'full');
    
    % factorize V using nmf_MDC_simple method

%% mdc test
    random_ = 1;
    if random_
        HI = abs(randn(size(HI)));
        WI = abs(randn(size(WI)));
        [ WMdcAscl1_2, HMdcAscl1_2, HRecord, E] = ...
            hyperNmfMdcAscl1_2(...
                V, HI, WI, ...
                0.001, ...
                20000, ...
                0.001, ...
                10 );
        WMdcRandomInit = WMdcAscl1_2;
        HMdcRandomInit = HMdcAscl1_2;
    else
        [ WMdcAscl1_2, HMdcAscl1_2, HRecord, E] = ...
            hyperNmfMdcAscl1_2(...
                V, HI, WI, ...
                0.001,... % tolObj
                20000, ... % maxIter
                0.001, ... % dDelta
                20 ... % fDelta
            );
        WMdcWellInit = WMdcAscl1_2;
        HMdcWellInit = HMdcAscl1_2;
    end   
    %% visualize rest
    VNmf = WMdcAscl1_2 * HMdcAscl1_2;
    figure;
    hold on
    bandIndx1 = 1;
    bandIndx2 = 2;
    scatter(V(:,bandIndx1), V(:,bandIndx2), 'c' ); hold on ; 
    scatter(HTrue(:, bandIndx1), HTrue(:, bandIndx2), 'filled', 'r');
    scatter(HI(:, bandIndx1), HI(:, bandIndx2), 'filled', 'g');
    scatter( VNmf(:,bandIndx1), VNmf(:,bandIndx2), 5, 'k' );
    scatter( HMdcAscl1_2(:,bandIndx1), HMdcAscl1_2(:,bandIndx2) , 'filled','k')
    plot(HRecord(:, 1, bandIndx1), HRecord(:, 1, bandIndx2), 'r-', 'MarkerSize', 5);
    plot(HRecord(:, 2, bandIndx1), HRecord(:, 2, bandIndx2), 'g-', 'MarkerSize', 5);
    plot(HRecord(:, 3, bandIndx1), HRecord(:, 3, bandIndx2), 'b-',  'MarkerSize', 5);
    plot(HRecord(:, 4, bandIndx1), HRecord(:, 4, bandIndx2), 'm-',  'MarkerSize', 5);