%% nmf_MDC_test with more than 3 endmembers
clear all
plot_ = false;
random_ = true;
expTimes = 50;
dataSize = 900;
noiseLevel = 0.001;
bandNum = 5;
emNum = 6;
emTrueData = zeros(emNum, bandNum);
emInitData = zeros(emNum, bandNum);
emMdcResult = zeros(emNum, bandNum);
emMvcResult = zeros(emNum, bandNum);
abunTrueData = zeros(dataSize, emNum);
abunInitData = zeros(dataSize, emNum);
abunMdcResult = zeros(dataSize, emNum);
abunMvcResult = zeros(dataSize, emNum);
sadMdc = zeros(emNum,expTimes);
sadMvc = zeros(emNum,expTimes);
for expIndex = 1:expTimes
    % generate true W, H, V
    emTrueData = abs( randn( emNum, bandNum ) );
    [hyperData, abunTrueData] = create4(dataSize, emTrueData);
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
    if random_==false
        emInitDataIndex = n_findr(hyperData, emNum);
        emInitData = hyperData(emInitDataIndex, :);

        % find initian W abundance using nmf
        % update abundance matrix only
        alpha = 1;
        tol = 0.1;
        maxIter = 5000;

        [abunInitData, E_I] = nmf_abundance(hyperData, emNum, emInitData,...
                            alpha, tol, maxIter);
        V_test = abunInitData * emInitData;
    end
    
    if plot_ == true
        figure(1)
        scatter(V_test(:,1), V_test(:,2), 'c', 'full'); hold on
        scatter(hyperData(:,1), hyperData(:,2), 'r')
        scatter(emInitData(:, 1), emInitData(:, 2), 'full');
    end
    
    % factorize V using nmf_MDC_simple method

%% mdc test
    if random_
        emInitData = abs(randn(emNum, bandNum));
        abunInitData = abs(randn(dataSize, emNum));
        [ abunMdcResult, emMdcResult, H_r, E] = ...
            hyperNmfMDC(...
                hyperData, emNum, abunInitData, emInitData, ...
                0.001, 0.3,...
                0.01, 30000 );
        W_test_random_init = abunMdcResult;
        H_test_random_init = emMdcResult;
    else
        [ abunMdcResult, emMdcResult, H_r, E] = ...
            hyperNmfMDC(...
                hyperData, emNum, abunInitData, emInitData, ...
                0.001, 0.3,...
                0.01, 30000 );
        W_test_well_init = abunMdcResult;
        H_test_well_init = emMdcResult;
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
            cur_sad = sad(emTrueData(em_j,:)', emMdcResult(em_i,:)');
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
        V_test = abunMdcResult * emMdcResult;
        figure(2)
        subplot(321)
        em_index1 = 1;
        em_index2 = 2;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(322)
        em_index1 = 1;
        em_index2 = 3;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(323)
        em_index1 = 1;
        em_index2 = 4;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);
    %     axis equal

         subplot(324)
        em_index1 = 1;
        em_index2 = 5;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);

         subplot(325)
        em_index1 = 2;
        em_index2 = 3;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);

         subplot(326)
        em_index1 = 2;
        em_index2 = 4;
        scatter(hyperData(:,em_index1), hyperData(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( V_test(:,em_index1), V_test(:,em_index2), 5, 'k' );
    %     scatter( H_true(:,1), H_true(:,2) , 100, 'filled','c')
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','k')
        plot(H_r(1, :,em_index1), H_r(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(H_r(2, :,em_index1), H_r(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(H_r(3, :,em_index1), H_r(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(H_r(4, :,em_index1), H_r(4, :,em_index2), 'm-',  'MarkerSize', 5);

        figure(3)
        subplot(121); plot(E_I)
        subplot(122); plot(E)
    end
    
 %% mvcnmf ref test
    [UU, SS, WW] = svd(hyperData');
    prin_comp = pca(hyperData);
    mean_data = mean(emTrueData, 1);

    [emMvcResult, abunMvcResult] = hyperNmfMVC(hyperData', emInitData', abunInitData', ... 
                            emTrueData', UU, prin_comp, mean_data, ...
                            0.001, 0.001, 30000, ...
                            0, 2, 1);
    V_mvc = emMvcResult * abunMvcResult;
    for em_i = 1:emNum
        tmp_sad = inf;
        for em_j = 1:emNum
            cur_sad = sad(emTrueData(em_j,:)', emMvcResult(:,em_i));
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
        scatter(emTrueData(:, em_index1), emTrueData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')

        subplot(232)
        em_index1 = 1;
        em_index2 = 3;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInitData(:, em_index1), emInitData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')

        subplot(233)
        em_index1 = 1;
        em_index2 = 4;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInitData(:, em_index1), emInitData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')

        subplot(234)
        em_index1 = 1;
        em_index2 = 5;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInitData(:, em_index1), emInitData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')

        subplot(235)
        em_index1 = 2;
        em_index2 = 3;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInitData(:, em_index1), emInitData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')

        subplot(236)
        em_index1 = 2;
        em_index2 = 4;
        scatter(hyperData(:, em_index1), hyperData(:, em_index2), 'c'); hold on;

        scatter(V_mvc(em_index1,:), V_mvc(em_index2, :), 5, 'k')
        scatter(emInitData(:, em_index1), emInitData(:, em_index2), 100, 'filled', 'm');
        scatter(emMvcResult(em_index1, :), emMvcResult(em_index2, :), 100, 'filled', 'k');
        scatter( emMdcResult(:,em_index1), emMdcResult(:,em_index2) , 100, 'filled','r')
    end
end

%% 
sadMdcAv = sum(sadMdc, 2)/expTimes;
sadMvcAv = sum(sadMvc, 2)/expTimes;
%-- plot(sadMdc/expTimes, 'r--*', 'LineWidth', 2); hold on
%-- plot(sadMvc/expTimes, '--*', 'LineWidth', 2)
%-- xlabel('Experimets', 'FontSize', 15, 'FontWeight', 'bold')
%-- ylabel('SAD', 'FontSize', 15, 'FontWeight', 'bold')
%-- legend('SSD Constraint','Volume Constraint', 'Orientation', 'Vertical')
hBar = bar([sadMdcAv, sadMvcAv]);
set(hBar(1), 'FaceColor', [0.2, 0.2, 0.2])
set(hBar(2), 'FaceColor', [0.6, 0.6, 0.6])
xlabel('endmember', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('SAD', 'FontSize', 15, 'FontWeight', 'bold')
legend('MDC-NMF', 'MVC-NMF', 'Orientation', 'Vertical')