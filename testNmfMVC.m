%% nmf_MDC_test with more than 3 endmembers
clear all
expTimes = 1;
dataSize = 900;
noiseLevel = 0.0;
bandNum = 3;
endNum = 3;
HTrue = zeros(endNum, bandNum);
HI = zeros(endNum, bandNum);
HMdc = zeros(endNum, bandNum);
WTrue = zeros(dataSize, endNum);
WI = zeros(dataSize, endNum);
WMdc = zeros(dataSize, endNum);
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
    [UU, SS, WW] = svd(V');
    prinComp = pca(V);
    mean_data = mean(HTrue, 1);

    [HMvc, WMvc] = hyperNmfMVC(V', HI', WI', ... 
                            HI', UU, prinComp, mean_data, ...
                            0.015, 0.001, 30000, ...
                            0, 2, 1);
    VMvc = HMvc * WMvc;  
    %% visualize rest
    VNmf = HMvc * WMvc;
    figure;
    emIndx1 = 1;
    emIndx2 = 2;
    scatter(V(:,emIndx1), V(:,emIndx2), 'c' ); hold on ; 
    scatter(HTrue(:, emIndx1), HTrue(:, emIndx2), 'filled', 'r');
    scatter(HI(:, emIndx1), HI(:, emIndx2), 'filled', 'b');
    scatter( VNmf(emIndx1, :), VNmf(emIndx2, :), 5, 'k' );
    scatter( HMvc(emIndx1,:), HMvc(emIndx2,:) , 'filled','k')