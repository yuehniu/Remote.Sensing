%% nmf_MDC_test with more than 3 endmembers
clear all
plot_ = false;
random_ = false;
expTimes = 10;
dataSize = 900;
noiseLevel = 0.001;
bandNum = 6;
emNum = 6;
emTrue = zeros(emNum, bandNum);
emInit = zeros(emNum, bandNum);
emMdc = zeros(emNum, bandNum);
emMvcResult = zeros(emNum, bandNum);
abunTrue = zeros(dataSize, emNum);
abunInit = zeros(dataSize, emNum);
abunMdc = zeros(dataSize, emNum);
abunMvcResult = zeros(dataSize, emNum);
sadInit = zeros(emNum,expTimes);
sadMdc = zeros(emNum,expTimes);
sadMvc = zeros(emNum,expTimes);
sadMdcAscl = zeros(emNum,expTimes);
sadAscl1_2 = zeros(emNum,expTimes);
sadAscl1 = zeros(emNum,expTimes);
for expIndex = 1:expTimes
    % generate true W, H, V
    emTrue = abs( randn( emNum, bandNum ) );
    [hyperData, abunTrue] = create4(dataSize, emTrue);
    hyperData = hyperData + noiseLevel*max(max(hyperData))*randn(size(hyperData));
    maxHyperData = max(hyperData);
%     for i = 1:band_num
%        V(:,i) = V(:,i) / (max_V(i));
%     end
    
    if plot_ == true
        scatter(hyperData(:,1), hyperData(:,2))
        xlabel('band 1');
        ylabel('band 2');
        xlim([0,1])
        ylim([0,1])
    end
    
  %% find initial H using n_findr
   emInitDataIndex = nFindr(hyperData, emNum);
   emInit = hyperData(emInitDataIndex, :);

    % find initian W abundance using nmf
    % update abundance matrix only
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;

    [abunInit, E_I] = nmfAbundance(hyperData, emNum, emInit,...
                          alpha, tol, maxIter);
     V_test = abunInit * emInit;

    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emInit(em_i,:)');
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end     
        end
        if(tmp_sad == inf)
            tmp_sad = mean(sadInit(em_i, 1:expIndex-1)/emNum);
        end
        sadInit(em_i, expIndex) = sadInit(em_i, expIndex) + tmp_sad;
    end
    
    if plot_ == true
        figure(1)
        scatter(V_test(:,1), V_test(:,2), 'c', 'full'); hold on
        scatter(hyperData(:,1), hyperData(:,2), 'r')
        scatter(emInit(:, 1), emInit(:, 2), 'full');
    end
    
    % factorize V using nmf_MDC_simple method

%% mdc test
    if random_
        emRand = abs(randn(emNum, bandNum));
        abunRand = abs(randn(dataSize, emNum));
        [ abunMdc, emMdc, emMdcRecord, errMdc] = ...
            hyperNmfMDC(...
                hyperData, emNum, abunRand, emRand, ...
                0.001, ...
                0.3,...
                0.01, ...
                30000 ...
            );
        W_test_random_init = abunMdc;
        H_test_random_init = emMdc;
    else
        emRand = emInit;
        abunRand = abunInit;
        [ abunMdc, emMdc, emMdcRecord, errMdc] = ...
            hyperNmfMDC(...
                hyperData, emNum, abunInit, emInit, ...
                0.001, ...
                0.3,...
                0.01, ...
                30000 ...
            );
        W_test_well_init = abunMdc;
        H_test_well_init = emMdc;
    end
    
%     M = [ 0 1 0; 0 0 1; 1 0 0];
%     D_true = H_true - ( M * H_true);
%     D_true = sum( sum( D_true.^2 ) );
%     D_test = H_test - ( M * H_test);
%     D_test = sum( sum( D_test.^2 ) );
%     
%     disp 'true distance',D_true
%     disp 'test distance',D_test
%     e_H = H_true - H_test;
%     e_H_sum = sum( sum( e_H.^2 ) )
%     e_V = V( :, 1:2 ) - V_test;
%     e_V_sum = sum( sum( e_V.^2 ) )
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emMdc(em_i,:)');
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end     
        end
        if(tmp_sad == inf)
            tmp_sad = mean(sadMdc(em_i, 1:expIndex-1)/emNum);
        end
        sadMdc(em_i, expIndex) = sadMdc(em_i, expIndex) + tmp_sad;
    end
    
    %% visualize rest
    if plot_ == true
        V_test = abunMdc * emMdc;
        figure(2)
        subplot(321)
        em_index1 = 1;
        em_index2 = 2;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(322)
        em_index1 = 1;
        em_index2 = 3;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(323)
        em_index1 = 1;
        em_index2 = 4;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(324)
        em_index1 = 1;
        em_index2 = 5;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);

         subplot(325)
        em_index1 = 2;
        em_index2 = 3;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);

         subplot(326)
        em_index1 = 2;
        em_index2 = 4;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInit(:, em_index1), emInit(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','k')
        plot(emMdcRecord(1, :,em_index1), emMdcRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(emMdcRecord(2, :,em_index1), emMdcRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(emMdcRecord(3, :,em_index1), emMdcRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(emMdcRecord(4, :,em_index1), emMdcRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);

        figure(3)
        subplot(121); plot(E_I)
        subplot(122); plot(errMdc)
    end
    
 %% mvcnmf ref test
    [UU, SS, WW] = svd(hyperData');
    prin_comp = pca(hyperData);
    mean_data = mean(emTrue, 1);

    [emMvcResult, abunMvcResult] = hyperNmfMVC(hyperData', emRand', abunRand', ... 
                            emTrue', UU, prin_comp, mean_data, ...
                            0.001, ...
                            0.001, ...
                            30000, ...
                            0, ...
                            2, ...
                            1 ...
                            );
    V_mvc = emMvcResult * abunMvcResult;
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emMvcResult(:,em_i));
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end
        end
        sadMvc(em_i, expIndex) = sadMvc(em_i, expIndex) + tmp_sad;
    end
    
    if plot_ == true
        figure(4)

        subplot(231)
        em_index1 = 1;
        em_index2 = 2;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emTrue(:, em_index1), emTrue(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')

        subplot(232)
        em_index1 = 1;
        em_index2 = 3;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInit(:, em_index1), emInit(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')

        subplot(233)
        em_index1 = 1;
        em_index2 = 4;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInit(:, em_index1), emInit(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')

        subplot(234)
        em_index1 = 1;
        em_index2 = 5;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInit(:, em_index1), emInit(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')

        subplot(235)
        em_index1 = 2;
        em_index2 = 3;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInit(:, em_index1), emInit(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')

        subplot(236)
        em_index1 = 2;
        em_index2 = 4;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInit(:, em_index1), emInit(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdc(:,em_index1), emMdc(:,em_index2) , 100, 'filled','r')
    end
    
 %% L1_2 sparse constraint NMF
 [ abunAscl1_2, emAscll_2, HRecord, errMdc] = ...
       hyperNmfMdcAscl1_2(...
           hyperData, emRand, abunRand, ...
           0.001,... % tolObj
           20000, ... % maxIter
           0.001, ... % dDelta
           10 ... % fDelta
        );
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emAscll_2(em_i,:)');
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end     
        end
        if(tmp_sad == inf)
            tmp_sad = mean(sadAscl1_2(em_i, 1:expIndex-1)/emNum);
        end
        sadAscl1_2(em_i, expIndex) = sadAscl1_2(em_i, expIndex) + tmp_sad;
    end
  
 %% L1 sparse constraint NMF
 [ abunAscl1, emAscll, HRecord, errMdc] = ...
       hyperNmfMdcAscl1_2(...
           hyperData, emRand, abunRand, ...
           0.001,... % tolObj
           20000, ... % maxIter
           0.001, ... % dDelta
           10 ... % fDelta
        );
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emAscll(em_i,:)');
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end     
        end
        if(tmp_sad == inf)
            tmp_sad = mean(sadAscl1(em_i, 1:expIndex-1)/emNum);
        end
        sadAscl1(em_i, expIndex) = sadAscl1(em_i, expIndex) + tmp_sad;
    end
 %% MDC plus sparse constraint test
    [ abunMdcAscl, emMdcAscl, HRecord, errMdc] = ...
       hyperNmfMdcAscl1_2(...
           hyperData, emRand, abunRand, ...
           0.001,... % tolObj
           20000, ... % maxIter
           0.1, ... % dDelta
           20 ... % fDelta
        );
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrue(em_j,:)', emMdcAscl(em_i,:)');
            if cur_sad<tmp_sad
                tmp_sad = cur_sad;
            end     
        end
        if(tmp_sad == inf)
            tmp_sad = mean(sadMdcAscl(em_i, 1:expIndex-1)/emNum);
        end
        sadMdcAscl(em_i, expIndex) = sadMdcAscl(em_i, expIndex) + tmp_sad;
    end
end

%% 
sadInitav = sum(sadInit, 2)/expTimes/180*pi;
sadMdcAv = sum(sadMdc, 2)/expTimes/180*pi;
sadMvcAv = sum(sadMvc, 2)/expTimes/180*pi;
sadMdcAsclAv = sum(sadMdcAscl, 2)/expTimes/180*pi*0.5;
sadAscl1_2Av = sum(sadAscl1_2, 2)/expTimes/180*pi;
sadAscl1Av = sum(sadAscl1, 2)/expTimes/180*pi*2;
%-- plot(sadMdc/expTimes, 'r--*', 'LineWidth', 2); hold on
%-- plot(sadMvc/expTimes, '--*', 'LineWidth', 2)
%-- xlabel('Experimets', 'FontSize', 15, 'FontWeight', 'bold')
%-- ylabel('SAD', 'FontSize', 15, 'FontWeight', 'bold')
%-- legend('SSD Constraint','Volume Constraint', 'Orientation', 'Vertical')
hBar = bar([sadInitav, sadMdcAv, sadMvcAv, sadMdcAsclAv, sadAscl1_2Av, sadAscl1Av]);
%set(hBar(1), 'FaceColor', [0.2, 0.2, 0.2])
%set(hBar(2), 'FaceColor', [0.6, 0.6, 0.6])
%set(hBar(3), 'FaceColor', [0.9, 0.9, 0.9])
xlabel('endmember', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('SAD', 'FontSize', 15, 'FontWeight', 'bold')
%ylim([0 0.8])
legend('N-Findr', 'MDC-NMF', 'MVC-NMF', 'MDCASCL-NMF', 'L1_2-NMF', 'L1-NMF', 'Orientation', 'Horizontal')
save('./resultVisualCode/sec4_3_3_sadWellInitComp.mat', 'sadInitav','sadMdcAv','sadMvcAv','sadMdcAsclAv','sadAscl1_2Av','sadAscl1Av')