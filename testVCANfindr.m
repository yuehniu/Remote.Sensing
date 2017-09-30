%% compare between MDC and MVC algorithm
    clear all
    % generate true W, H, V
    sampleNum = 900;
    noiseLevel = 0.2;
    bandNum = 5;
    endNum = 3;
    
    % compare several times
    HTrue = abs( randn( endNum, bandNum ) );
    V = create4(sampleNum, HTrue);
    
    % run VCA
    HVca = hyperVca(V', endNum);
    
    % run Nfindr
    indxNfindr = nFindr(V, endNum);
    HNfindr = V(indxNfindr, :);
    
    % visualize result
    figure;
    hold on
    scatter(V(:,1), V(:,2));
    scatter(HTrue(:, 1), HTrue(:, 2), 'filled', 'r');
    scatter(HVca(1,:), HVca(2,:), 80, 'filled', 'k')
    scatter(HNfindr(:, 1), HNfindr(:, 2), 'filled', 'm')