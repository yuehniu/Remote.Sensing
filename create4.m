function [R, W]=create4(M,E)
% M: number of hyper data
% E: endmember data
    [N band]=size(E);%E=rand(N,band);
    W = zeros(M, N);
    for i=1:M
        l(1:N)=1;
        a=drchrnd(l,1);
        while max(a)>0.8
            a=drchrnd(l,1);
        end
        R(i,:)=a*E;
        W(i,:) = a;
    end
end

function r = drchrnd(a,n)
    % take a sample from a dirichlet distribution
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
end