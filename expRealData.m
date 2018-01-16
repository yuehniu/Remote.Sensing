%% expriment on real data.
    clear all
%% read real hyper data
    plot_ = false;
    
    hyperData = imread('./data/cup95eff.tif');

    % rearange hyper data
    hyperDataLine = permute(hyperData, [3, 1, 2]);
    sizeData = size(hyperDataLine);
    hyperDataLine = reshape(hyperDataLine,...
    sizeData(1), sizeData(2)*sizeData(3));
    hyperDataLine = double(hyperDataLine');
    maxOrigData = max(max(hyperDataLine));
    hyperDataLine = hyperDataLine / max(max(hyperDataLine));
    
    % pca analysis
    bandNum = 10;
    emNum = 10;
    scal = 0.1;
    data_size = sizeData(2) * sizeData(3);
    [coeff,score,~,~,explained,mu] = pca(...
      hyperDataLine', ...
      'NumComponents', bandNum,...
      'Centered', false);
    minCoeff = min(coeff);
    maxCoeff = max(coeff);
    mainDataLine = coeff - ...
        repmat(minCoeff, data_size, 1) + ...
        repmat(scal*maxCoeff, data_size, 1); 
    maxMainData = max(mainDataLine);
    mainDataLine = mainDataLine / diag(maxMainData);
       
    %% find initial H using n_findr
    tic;
    emInitIndex = nFindr(mainDataLine, emNum);
    toc;
    emInitData = mainDataLine(emInitIndex, :);
    
    % find initian W abundance using nmf
    % update abundance matrix only
    alpha = 1;
    tol = 0.1;
    maxIter = 5000;
    tic;
    [abunInitData, E_I] = nmfAbundance(mainDataLine, emNum, emInitData,...
                        alpha, tol, maxIter);
    toc;
    VNfindr = abunInitData * emInitData;
    
    %% visualize init endmember
    if plot_
        figure
        show_index1 = 2;
        show_index2 = 5;
        point_stride = 1;
        scatter(VNfindr(1:point_stride:data_size,show_index1),...
                VNfindr(1:point_stride:data_size,show_index2), ...
                'c', 'full'); hold on
        scatter(mainDataLine(1:point_stride:data_size,show_index1),...
                mainDataLine(1:point_stride:data_size,show_index2), ...
                'r')
        scatter(emInitData(:, show_index1), emInitData(:, show_index2), 'k', 'full');
    end
    
    %% mdc test
    tic
%     H_I = abs(randn(size(H_I)));
%     W_I = abs(randn(size(W_I)));
    [ WMdc, HMdc, HRecord, E] = ...
        hyperNmfMDC(...
            mainDataLine, emNum, abunInitData, emInitData, ...
            0.5, 0.01,...
            90, 30000 );
     toc   
    %% visualize rest
    if plot_
        VNfindr = WMdc * HMdc;
        figure
        em_index1 = 2;
        em_index2 = 3;
        scatter(mainDataLine(:,em_index1), ...
                mainDataLine(:,em_index2), 'c' ); hold on ; 
        scatter(emInitData(:, em_index1), emInitData(:, em_index2));
        scatter( VNfindr(:,em_index1), VNfindr(:,em_index2), 'k' );
        scatter( HMdc(:,em_index1), HMdc(:,em_index2) , 100, 'filled','r')
        plot(HRecord(1, :,em_index1), HRecord(1, :,em_index2), 'r-', 'MarkerSize', 5);
        plot(HRecord(2, :,em_index1), HRecord(2, :,em_index2), 'g-', 'MarkerSize', 5);
        plot(HRecord(3, :,em_index1), HRecord(3, :,em_index2), 'b-',  'MarkerSize', 5);
        plot(HRecord(4, :,em_index1), HRecord(4, :,em_index2), 'm-',  'MarkerSize', 5);
    end
    
    %% MDSC test
     [ WMdcAscl1_2, HMdcAscl1_2, HRecord, E] = ...
       hyperNmfMdcAscl1_2(...
          mainDataLine, emInitData, abunInitData, ...
          68,... % tolObj
          20000, ... % maxIter
          0.001, ... % dDelta
          20 ... % fDelta
      );
    
    %% compare H result and H init
    % N-Findr
    HNfindrOrig = emInitData * diag(maxMainData);
    HNfindrOrig = HNfindrOrig + ...
        repmat(minCoeff, emNum, 1) - ...
        repmat(scal*maxCoeff, emNum, 1);
    HNfindrOrig = score * HNfindrOrig';
    %H_I_orig = H_I_orig + repmat(mu(H_init_index), 50, 1);
    HNfindrOrig = HNfindrOrig * maxOrigData;
    
    % MDC
    HMdcOrig = HMdc * diag(maxMainData);
    HMdcOrig = HMdcOrig + ...
        repmat(minCoeff, emNum, 1) - ...
        repmat(scal*maxCoeff, emNum, 1);
    HMdcOrig = score * HMdcOrig';
    %H_test_orig = H_test_orig + repmat(mu(H_init_index), 50, 1);
    HMdcOrig = HMdcOrig * maxOrigData;
    
    % MDSC
    HMdcAscl1_2Orig = HMdcAscl1_2 * diag(maxMainData);
    HMdcAscl1_2Orig = HMdcAscl1_2Orig + ...
        repmat(minCoeff, emNum, 1) - ...
        repmat(scal*maxCoeff, emNum, 1);
    HMdcAscl1_2Orig = score * HMdcAscl1_2Orig';
    %H_test_orig = H_test_orig + repmat(mu(H_init_index), 50, 1);
    HMdcAscl1_2Orig = HMdcAscl1_2Orig * maxOrigData;
    
    %% visualize abundance
    for i = 1:emNum
        cur_abundance = WMdc(:, i);
        cur_abundance = reshape(cur_abundance, sizeData(2), sizeData(3));
        figure(i)
        imshow(cur_abundance, 'Border', 'tight')
    end
    
%% ref mvc test
%     [UU, SS, WW] = svd(hyper_data_line');
%     prin_comp = pca(hyper_data_line);
%     H_I = abs(randn(size(H_I)));
%     W_I = abs(randn(size(W_I)));
    mean_data = mean(emInitData, 1);

    [H_mvc, W_mvc] = mvcnmf(mainDataLine', emInitData', abunInitData', ... 
                            emInitData', score', mainDataLine', mean_data, ...
                            0.015, 35, 50, ...
                            0, 2, 1);
    V_mvc = H_mvc * W_mvc;
    
 %% visualize mvc test
    V_mvc = W_mvc' * H_mvc';
    figure(2)
    em_index1 = 2;
    em_index2 = 6;
    scatter(mainDataLine(:,em_index1), ...
            mainDataLine(:,em_index2), 'c' ); hold on ; 
    scatter(emInitData(:, em_index1), emInitData(:, em_index2));
    scatter( V_mvc(:,em_index1), V_mvc(:,em_index2), 'k' );
    scatter( H_mvc(em_index1, :), H_mvc(em_index2, :) , 100, 'filled','r')