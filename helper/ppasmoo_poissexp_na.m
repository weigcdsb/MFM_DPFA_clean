function [x,W,lam,lamPred] = ppasmoo_poissexp_na(n,C,d,x0,W0,A,b,Q)

% to debug
% n = Y_train;
% C = C_fit(:,:,1);
% d = d_tmp;
% x0 = x0_fit(:,1);
% W0 = Q0;
% A = A_fit(:,:,1);
% b = b_fit(:,1);
% Q = Q_fit(:,:,1);


warning('off');
lastwarn('')
nCell = size(n, 1);
T = size(n, 2);

% Preallocate
x   = zeros(length(x0),T);
W   = zeros([size(W0) T]);
lamPred = zeros(nCell, T);

% Initialize
x(:,1)   = x0;
W(:,:,1) = W0;
lamPred(:,1)   = exp(C*x0 + d);

xpred = x;
Wpred = W;
lam = lamPred;

I = eye(size(W0));

% Forward-Pass (Filtering)
for i=2:size(n,2)
    xpred(:,i) = A*x(:,i-1) + b;
    lamPred(:,i) = exp(C*xpred(:,i) + d);
    Wpred(:,:,i) = A*W(:,:,i-1)*A' + Q;
    
    INFO = zeros(size(W0));
    SCORE = zeros(size(x0));
    
    for k=1:nCell
        if(~isnan(n(k,i)))
            INFO = INFO + C(k,:)'*(lamPred(k,i))*C(k,:);
            SCORE = SCORE + C(k,:)'*(n(k,i) - lamPred(k, i));
        end
    end
    Wpostinv = inv(Wpred(:,:,i)) + INFO;
    W(:,:,i) = inv(Wpostinv);
    W(:,:,i) = (W(:,:,i) + W(:,:,i)')/2;
    
    x(:,i)  = xpred(:,i) + W(:,:,i)*SCORE;
    
    lam(:,i) = exp(C*x(:,i) + d);
    
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        error('Error: singular matrix');
    end
end
% lastwarn('')


% Backward-Pass (RTS)
for i=(T-1):-1:1
    Wi = inv(Wpred(:,:,i+1));
    J = W(:,:,i)*A'*Wi;
    x(:,i) = x(:,i) + J*(x(:,i+1) - xpred(:,i+1));
    W(:,:,i) = W(:,:,i) + J*(W(:,:,i+1)-Wpred(:,:,i+1))*J';
    
    W(:,:,i) = (W(:,:,i) + W(:,:,i)')/2;
    lam(:,i) = exp(C*x(:,i) + d);
end
warning('on');


end