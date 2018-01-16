%% nmf_MDC_test with more than 3 endmembers
    clear all
exp_times = 50;
data_size = 900;
noise_level_arr = [0.3162, 0.1778, 0.1, 0.0566, 0.0316];
noise_level_arr = [0];
band_num = 5;
N = 6;
n_N = 5;
sad_em_result_exp = zeros(3, N);
sad_ab_result_exp = zeros(3, N);
sad_em_with_noise = zeros(3, n_N);
sad_ab_with_noise = zeros(3, n_N);

for noise_i = 1:1
    noise_level = noise_level_arr(noise_i);
for exp_i = 1:exp_times
    % generate true W, H, V
    W_true = abs( randn( data_size, N));
    sum_W = sum( W_true, band_num );
    for i = 1:data_size;
        W_true( i, : ) = W_true( i,: ) ./ sum_W(i);
    end
    H_true = abs( randn( N, band_num ) );
    V = W_true * H_true;  
    [V, W_true] = create4(data_size, H_true);
    V = V + noise_level*rand(data_size, band_num);
    
    % find initial H using n_findr
    disp([num2str(exp_i), 'th iteration: init..'])
    H_init_index = nFindr(V, N);
    H_I = V(H_init_index, :);
    
    % find initian W abundance using nmf
    % update abundance matrix only
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;
    [W_I, E_I] = nmfAbundance(V, N, H_I,...
                        alpha, tol, maxIter);
    V_test = W_I * H_I;
    
    % factorize V using nmf_MDC_simple method

%% mdc test
    random_ = 0;
    disp([num2str(exp_i), 'th iteration: NMF_MDC..'])
    if random_
        H_I = abs(randn(size(H_I)));
        W_I = abs(randn(size(W_I)));
        [ W_test, H_test, H_r, E] = ...
            hyperNmfMDC(...
                V, N, W_I, H_I, ...
                0.001, 0.3,...
                0.01, 30000 );
        W_test_random_init = W_test;
        H_test_random_init = H_test;
    else
        [ W_test, H_test, H_r, E] = ...
            hyperNmfMDC(...
                V, N, W_I, H_I, ...
                0.001, 0.3,...
                0.01, 30000 );
        W_test_well_init = W_test;
        H_test_well_init = H_test;
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
    
    % visualize rest
    V_test = W_test * H_test;
    
 %% mvcnmf ref test
    [UU, SS, WW] = svd(V');
    prin_comp = pca(V);
    mean_data = mean(H_true, 1);
    disp([num2str(exp_i), 'th iteration: NMF_MVC..'])
%     W_I_mvc = abs(randn(size(W_I)));
    [H_mvc, W_mvc] = hyperNmfMVC(V', H_I', W_I', ... 
                            H_true', UU, prin_comp, mean_data, ...
                            0.015, 0.001, 30000, ...
                            0, 2, 1);
    V_mvc = H_mvc * W_mvc;

    
    %% calculate sad between vectors
    em_num = N;
    sad_em_result = zeros(3, em_num);
    sad_ab_result = zeros(3, em_num);
    em_I = zeros(1,em_num);
    em_indx_mdc = zeros(1,em_num);
    em_indx_mvc = zeros(1,em_num);
    % find the closed real em
    tmp_sad = inf;
    for i = 1:em_num
        tmp_sad = inf;
        for j = 1:em_num
            if(sad(H_true(i,:)', H_I(j,:)') < tmp_sad)
                em_I(i) = j;
                tmp_sad = sad(H_true(i,:)', H_I(j,:)');
            end
        end
    end
    tmp_sad = inf;
    for i = 1:em_num
        tmp_sad = inf;
        for j = 1:em_num
            if(sad(H_true(i,:)', H_test(j,:)') < tmp_sad)
                em_indx_mdc(i) = j;
                tmp_sad = sad(H_true(i,:)', H_test(j,:)');
            end
        end
    end

    for i = 1:em_num
        tmp_sad = inf;
        for j = 1:em_num
            if(sad(H_true(i,:)', H_mvc(:,j)) < tmp_sad)
                em_indx_mvc(i) = j;
                tmp_sad = sad(H_true(i,:)', H_mvc(:,j));
            end
        end
    end

    compare_term = 1; % 1: endmember 0: abundance

    for i = 1:em_num
        sad_em_result(1, i) = sad(H_I(em_I(i),:)', H_true(i,:)');
    end

    for i = 1:em_num
        sad_em_result(2, i) = sad(H_test(em_indx_mdc(i),:)', H_true(i,:)');
    end

    for i = 1:em_num
        sad_em_result(3, i) = sad(H_mvc(:, em_indx_mvc(i)), H_true(i, :)');
    end
    for i = 1:em_num
        sad_ab_result(1, i) = sad(W_I(:, em_I(i)), W_true(:, i));
    end

    for i = 1:em_num
        sad_ab_result(2, i) = sad(W_test(:,em_indx_mdc(i)), W_true(:, i));
    end

    for i = 1:em_num
        sad_ab_result(3, i) = sad(W_mvc(em_indx_mvc(i), :)', W_true(:, i));
    end
    sad_em_result_exp = sad_em_result_exp + sad_em_result;
    sad_ab_result_exp = sad_ab_result_exp + sad_ab_result;
end

sad_em_result_exp = sad_em_result_exp / exp_times;
sad_ab_result_exp = sad_ab_result_exp / exp_times;

sad_em_with_noise(:, noise_i) = mean(sad_em_result_exp, 2);
sad_ab_with_noise(:, noise_i) = mean(sad_ab_result_exp, 2);
end