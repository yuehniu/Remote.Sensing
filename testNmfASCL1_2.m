%% hyperNmfASCL1_2 test with more than 3 endmembers
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 3;
endNum = 3;
HTrue = zeros(endNum, bandNum);
HI = zeros(endNum, bandNum);
HASCL1_2 = zeros(endNum, bandNum);
WTrue = zeros(dataSize, endNum);
WI = zeros(dataSize, endNum);
WASCL1_2 = zeros(dataSize, endNum);
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
%     xlim([0,1])
%     ylim([0,1])
    
    %% find initial H using n_findr
    HInitIndx = nFindr(V, endNum);
    HI = V(HInitIndx, :);
    
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;
    [WI, EI] = nmfAbundance(V, endNum, HI,...
                        alpha, tol, maxIter);
    VASCL1_2 = WI * HI;
    
    figure;
    scatter(VASCL1_2(:,1), VASCL1_2(:,2), 'c', 'full'); hold on
    scatter(V(:,1), V(:,2), 'r')
    scatter(HI(:, 1), HI(:, 2), 'full');
    
    % factorize V using nmf_MDC_simple method

%% ASCL1_2 test
    random_ = 0;
    if random_
        HI = abs(randn(size(HI)));
        WI = abs(randn(size(WI)));
        [ WASCL1_2, HASCL1_2, HRecord, E] = ...
            hyperNmfMDC(...
                V, endNum, WI, HI, ...
                0.001, 0.3,...
                0.01, 30000 );
        WMdcRandomInit = WASCL1_2;
        HMdcRandomInit = HASCL1_2;
    else
        [ WASCL1_2, HASCL1_2, HRc] = ...
            hyperNmfASCL1_2(...
                V', HI', WI',...
                0.0001,... % tolObj
                300,... % maxIter
                0.05 ... %fDelta
                );
        WMdcWellInit = WASCL1_2;
        HMdcWellInit = HASCL1_2;
    end   
    %% visualize rest
    VASCL1_2 = HASCL1_2 * WASCL1_2;
    figure;
    bandIndx1 = 1;
    bandIndx2 = 2;
    scatter(V(:,bandIndx1), V(:,bandIndx2), 'c' ); hold on ; 
    scatter(HTrue(:, bandIndx1), HTrue(:, bandIndx2), 'filled', 'r');
    scatter(HI(:, bandIndx1), HI(:, bandIndx2), 'filled', 'b');
    scatter( VASCL1_2(bandIndx1,:), VASCL1_2(bandIndx2,:), 5, 'k' );
    scatter( HASCL1_2(bandIndx1, :), HASCL1_2(bandIndx2,:) , 'filled','k')
    plot(HRc(1, :,bandIndx1), HRc(1, :,bandIndx2), 'r-', 'MarkerSize', 5);
    plot(HRc(2, :,bandIndx1), HRc(2, :,bandIndx2), 'g-', 'MarkerSize', 5);
    plot(HRc(3, :,bandIndx1), HRc(3, :,bandIndx2), 'b-',  'MarkerSize', 5);