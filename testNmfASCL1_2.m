%% hyperNmfASCL1_2 test with more than 3 endmembers
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 4;
emNum = 4;
HTrue = zeros(emNum, bandNum);
HI = zeros(emNum, bandNum);
HASCL1_2 = zeros(emNum, bandNum);
WTrue = zeros(dataSize, emNum);
WI = zeros(dataSize, emNum);
WASCL1_2 = zeros(dataSize, emNum);
% for i = 1:exp_times
    % generate true W, H, V
    HTrue = abs( randn( emNum, bandNum ) );
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
    HInitIndx = nFindr(V, emNum);
    HI = V(HInitIndx, :);
    
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;
    [WI, EI] = nmfAbundance(V, emNum, HI,...
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
        % HI = abs(randn(size(HI)));
        % WI = abs(randn(size(WI)));
        WI = WI ./ ( repmat( sum(WI,2),1, emNum ) );
        [ WASCL1_2, HASCL1_2, HRcL1_2, errRcL1_2, objRcL1_2] = ...
            hyperNmfASCL1_2(...
                V', HI', WI',...
                0.001,... % tolObj
                20000,... % maxIter
                10 ... %fDelta
                );
        WMdcRandomInit = WASCL1_2;
        HMdcRandomInit = HASCL1_2;
    else
        [ WASCL1_2, HASCL1_2, HRcL1_2, errRcL1_2, objRcL1_2] = ...
            hyperNmfASCL1_2(...
                V', HI', WI',...
                0.001,... % tolObj
                20000,... % maxIter
                20 ... %fDelta
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
    
    for i = 1:emNum
        plot(HRcL1_2(i, :,bandIndx1), HRcL1_2(i, :,bandIndx2), 'r-.', 'MarkerSize', 5);
    end
    
    figure;
    hold on;
    plot(errRcL1_2, 'r');
    plot(objRcL1_2, 'k')
    plot(objRcL1_2-errRcL1_2, 'm')