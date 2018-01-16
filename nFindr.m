%% N-FINDR function for endmember extraction
%-----------------------------------------------------------------------------------
    % this function is a realization for a new endmember extraction method, which is
    % similiar to the original N-FINDR method.
%-----------------------------------------------------------------------------------
function [endmember_index, V] = nFindr (hyspectral_data, endmember_number)
    % hyperspectral_data:   input hyperspectral_data
    % endmember_number:     input endmember endmember_number
    % endmember_index:      output index of endmember in original hyperspectral_data
    % V:                    output volumn conposed by endmember
    [row_hyper, col_hyper] = size (hyspectral_data);
    A = zeros (col_hyper, endmember_number); % store diff vector
    V_max = 0; 
    endmember_idx = zeros(1,endmember_number);
    endmember_vector = zeros (col_hyper, endmember_number); % store endmember_vector in every iteration
    
    flag_skip = 0;
    [man_val, min_indx]  = min( hyspectral_data(:,1));
    endmember_idx(1) = min_indx;
    endmember_vector(:,1) = hyspectral_data(min_indx, :)';
    
    disp_str = ['1th endmember found: [', num2str(endmember_vector(:,1)') ']'];
    disp(disp_str);
    
    for i = 2:endmember_number % iteration for endmember
        A = zeros (col_hyper,i-1);
        if i==1 % the first endmember extraction
            for j = 1:row_hyper % traverse all hyperspectral data
                A = hyspectral_data(j,:)';
                AAT = det (A' * A);
                AAT = sqrt (AAT); % volumn of endmember
                if AAT > V_max % update endmember
                    V_max = AAT;
                    endmember_vector(:,i) = A(:,i);
                    endmember_idx(i) = j;
                end
            end
        else
            for j = 1:row_hyper
                for jj = 1:i 
                    if (j==endmember_idx(jj)) || sum(hyspectral_data(j,:)'-endmember_vector(:,jj))==0
                        flag_skip = 1;
                        break;
                    end
                end
                if flag_skip==1
                    flag_skip = 0;
                    continue;
                end
                for ii = 1:i-1
                    A(:,ii) = hyspectral_data(j,:)' - endmember_vector(:,ii);
                end
                AAT = det (A' * A);
                AAT = sqrt (AAT) / factorial(i-1);
                if AAT > V_max
                    V_max = AAT;
                    endmember_vector(:,i) = hyspectral_data(j,:)';
                    endmember_idx(i) = j;
                end
                flag_skip = 0;
            end
        end
        V_max = 0;
        
        disp_str = [num2str(i) 'th endmember found: [', num2str(endmember_vector(:,i)') ']'];
        disp(disp_str);
    end
    endmember_index = endmember_idx;
    V = V_max;   
end
