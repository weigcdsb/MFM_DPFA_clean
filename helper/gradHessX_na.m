function gradHess = gradHessX_na(vecX, d_tmp, C_trans_tmp, x0_tmp,...
    Q0_tmp, Q_tmp, A_tmp, b_tmp, Y_train)

% to debug
% vecX = X_tmp(:);
% d_tmp = d_tmp;
% C_trans_tmp = C_fit(:,:,1);
% x0_tmp = x0_fit(:,1);
% Q0_tmp = Q0;
% Q_tmp = Q_fit(:,:,1);
% A_tmp = A_fit(:,:,1);
% b_tmp = b_fit(:,1);
% Y_train = Y_train;

T = size(Y_train,2);
X_all = reshape(vecX, [], T);
lam_tmp = exp(C_trans_tmp*X_all + d_tmp);

hessup = repmat((Q_tmp\A_tmp)', 1, 1, T-1);
hessub = repmat(Q_tmp\A_tmp, 1, 1, T-1);
hessmed = repmat(zeros(size(X_all, 1)),1,1,T);
score = zeros(size(X_all));
for t = 1:T
    
    obsIdx = find(~isnan(Y_train(:,t)));
    score(:,t) = C_trans_tmp(obsIdx,:)'*(Y_train(obsIdx,t) - lam_tmp(obsIdx,t));
    
    if (t==1)
        hessmed(:,:,t) = -C_trans_tmp(obsIdx,:)'*...
            diag(lam_tmp(obsIdx,t))*C_trans_tmp(obsIdx,:)-...
            inv(Q0_tmp)- A_tmp'*(Q_tmp\A_tmp);
    elseif (t == T)
        hessmed(:,:,t) = -C_trans_tmp(obsIdx,:)'*diag(lam_tmp(obsIdx,t))*...
            C_trans_tmp(obsIdx,:) -inv(Q_tmp);
    else
        hessmed(:,:,t) = -C_trans_tmp(obsIdx,:)'*diag(lam_tmp(obsIdx,t))*...
            C_trans_tmp(obsIdx,:) -inv(Q_tmp) - A_tmp'*(Q_tmp\A_tmp);
    end
end

gradHess{1} = score + [-Q0_tmp\(X_all(:,1) - x0_tmp)+...
    A_tmp'*(Q_tmp\(X_all(:,2) - A_tmp*X_all(:,1)-b_tmp)),...
    -Q_tmp\(X_all(:,2:(T-1)) - A_tmp*X_all(:,1:(T-2))-b_tmp)+...
    A_tmp'*(Q_tmp\(X_all(:,3:T) - A_tmp*X_all(:,2:(T-1))-b_tmp)),...
    -Q_tmp\(X_all(:,T) - A_tmp*X_all(:,T-1)-b_tmp)];
gradHess{1} = gradHess{1}(:);

gradHess{2} = blktridiag(hessmed,hessub,hessup);
% gradHess{2} = (gradHess{2} + gradHess{2}')/2;


end