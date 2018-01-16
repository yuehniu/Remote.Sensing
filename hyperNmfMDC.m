%% Minimum distant constraint for 
%  nonnegative matrix factorization
function [ W, H, H_r, E] = hyperNmfMDC(...
    V, N, W_I, H_I, ...
    lamda0, lamda1,... 
    Tolerance, IterMax)
% W,H:       output for factorized matrixs
% V:         original lined hyper data with some bands
% N:        endmember number
% alpha0:    distance constraint factor
% alpha1:    sum to 1 constraint factor
% Tolerance: max error when to stop iteration
% InterMax:  max interation number for stop iteration

%% start process
    M = eye(N);
    M(N, :) = [1, zeros(1, N-1)];
    P = M' + M;
    I = eye(N);
    [ row_V, col_V ] = size( V );
    V = [ V, lamda1*ones(row_V, 1) ];
    W = W_I;
%     sum_W = sum(W, 2);
%     for i = 1:row_V
%         W(i, :) = W(i, :) ./ sum_W(i);
%     end
    H = H_I;
    H = [ H lamda1*ones(N,1) ];
    
    % record iteration for H
    H_record = zeros(N, IterMax, col_V);
    % calculate linear least square error
    E = zeros(1, IterMax);
    e = (V - (W * H));
    e_sum = sum( sum( e.^2 ) );
    iter_cur = 1;
    E(1, iter_cur) = e_sum;
    for i = 1:N
        H_record(i, iter_cur, :) = H(i, 1:col_V);
    end
    % start iteration
    while( e_sum > Tolerance )
        if( iter_cur == IterMax )
            break;
        end
        
        % update H
        normi = (W' * V) + (lamda0 * P * H);
        denormi = ( (W' * W) + (2*lamda0 * I) ) * H;
        update_factor = normi ./ denormi;
        H = H .* update_factor;
        H( :,col_V+1 ) = lamda1 * ones( N, 1 );
        
        % update W
        normi = V * H';
        denormi = W * H * H';
        update_factor = normi ./ denormi;
        W = W .* update_factor;
        
        % calculate error
        e = (V - ( W * H ));
        e_sum = sum( sum( e.^2 ) );
        
        disp_str = [ '[' num2str(iter_cur), ']', ...
                    ' Loss: ' num2str(e_sum), ...
                    ' Initial Distance: ' num2str( norm(H_I - M*H_I, 'fro') ), ...
                    ' Estimate Distance: ' num2str( norm(H - M*H, 'fro') ) ];
        disp(disp_str);
        
        iter_cur = iter_cur + 1;
        E(1, iter_cur) = e_sum;
        
        for i = 1:N
            H_record(i, iter_cur, :) = H(i, 1:col_V);
        end
    end
    
    H = H( :,1:col_V );
    H_r = H_record(:, 1:iter_cur, :);
    E = E(1, 1:iter_cur);
end

