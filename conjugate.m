function S = conjugate(X,A,S,maxiter,U,meanData,tao)
%
% Projected conjugate gradient learning
%

% initialize
maxiter = 1;
AtA = A'*A;
AtX = A'*X;
[L,N] = size(X); [c,N] = size(S);

% the constraint is only included when estimating A
cons = 0;
if L>N
    cons = 1;
    meanData = meanData'*ones(1,c);
    C = [ones(1,c); zeros(c-1,c)];
    B = [zeros(1,c-1); eye(c-1)];
    Z = C+B*U'*(S'-meanData);
    ZD = pinv(Z)*B*U';
    detz2 = det(Z)^2;
end

% initial gradient
if cons == 1
    gradp = AtA*S - AtX + tao*detz2*ZD;
else
    gradp = AtA*S - AtX;
end
        
% initial conjugate direction
conjp = gradp;
S = S - 0.001*gradp;

for iter = 1:maxiter
    
    % new gradient
    if cons == 1
        Z = C+B*U'*(S'-meanData);
        ZD = pinv(Z)*B*U';
        detz2 = det(Z)^2;
        grad = AtA*S - AtX + tao*detz2*ZD;
    else
        grad = AtA*S - AtX;
    end
    
    % parameter beta
    beta = sum(grad.*(grad-gradp),1)./(sum(gradp.^2,1));
    
    % new conjugate direction
    conj = -grad+ones(c,1)*beta.*conjp;
    
    % steplength: line search to find exact minimizer of alpha
    AAd = AtA*conj;
    alpha = sum(conj.*(-grad),1)./max(sum(conj.*AAd,1),eps);
       
    % update
    S = S + conj.*repmat(alpha,size(AtA,2),1);   
    
    % save the current directions
    gradp = grad;
    conjp = conj;
    
end

