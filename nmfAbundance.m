function [W, E, IT] = nmfAbundance (V, N, H_I, alpha, tolerance, maxIter)
%% this function use nmf-similiar method to factorization hyperspectrum
% W,H:          output----factored output for hyperspectral data
% V:            input ----original hyperspectroal
% N:            input ----endmember numbers
% H_I:			H initial value
% alpha:        input ----scale factor 
% tolerance:    input ----input for error tolerance
% MaxIter:      input ----input for max iteration times

%%
[row_V, col_V] = size ( V );
V_EXT = [alpha*V, ones(row_V,1)];

V_cut = V_EXT(1:1:row_V,:);
[row_V, col_V] = size ( V_cut );
W = abs( randn ( row_V, N ));
sum_W = sum (W, 2);
for i = 1:row_V;
    W(i,:) = W(i,:) ./ sum_W(i);
end
% Hi = randi ( row_V );
H = [H_I,ones(N,1)]; % add sum is 1 to the equation
% for i = 1:N
%     Hi = randi ( row_V );
%     H(i,:) = V_cut(Hi,:);
% end

j = 1; % recorde iter

% record error for every iteration
E = zeros(1,maxIter);
% get the error

e = V_cut - W * H;
e2 = sum ( sum (e .^2) ) ;
Iter_n = 1; % iteration times

    while ( e2 > tolerance )
        if (Iter_n == maxIter )
            break;
        end
        
        
        % get nex W
        VH = V_cut * H';
        WHH = W * ( H * H' );
        W = W .* ( VH ./ WHH ); % update W
%         W_sum = sum(W, 2);
%         [r_, c_] = size(W);
%         for i_ = 1:c_
%             W(:,i_) = W(:,i_) ./ W_sum;
%         end

        % get error
        e = V_cut - W * H;
        e2 = sum ( sum (e .^2) );
        
        j = j + 1;       
        E(Iter_n) = e2;
            
        disp_str = ['[' num2str(Iter_n) ']',...
                    ' Loss: ' num2str(e2)];
        disp(disp_str);
        Iter_n = Iter_n + 1;
    end
   
    E = E(1:Iter_n);
    IT = Iter_n;
end