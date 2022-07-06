function [theta_b, C_b] = bdmcmc(theta_a, C_a, delt_tmp, birth_rate,...
    birth_time,alpha, Y_tmp, lab_unique, lab_all, prior, p_max)

% to debug
% theta_a = THETA{g};
% C_a = C_fit(:,:,g);
% delt_tmp = delt_fit(:,g);
% birth_rate = birth_rate;
% birth_time = birth_time;
% alpha = alpha;
% Y_tmp = Y;
% nClus = nClus;
% lab_tmp = Lab;


T = size(Y_tmp,2);
theta_b = theta_a;
C_b = C_a;

for j = lab_unique(:)'
    t_fa = 0;
    obsIdx = find(lab_all == j);
    N_tmp = length(obsIdx);
    
    while(t_fa < birth_time)
        
        mll_del = zeros(theta_b(j).p,1);
        for kk = 1:theta_b(j).p
            colSelect = setdiff(1:theta_b(j).p, kk);
            mll_del(kk) = poiLogMarg(reshape(Y_tmp(obsIdx,:)', [], 1),...
                repmat(theta_b(j).X(colSelect,:)', N_tmp,1),...
                kron(delt_tmp(obsIdx), ones(T,1)) +...
                repmat(theta_b(j).mu', N_tmp,1));
        end
        
        mll_all = poiLogMarg(reshape(Y_tmp(obsIdx,:)', [], 1),...
            repmat(theta_b(j).X', N_tmp,1),...
            kron(delt_tmp(obsIdx), ones(T,1)) +...
            repmat(theta_b(j).mu', N_tmp,1));
        
        delta_j = exp(mll_del + log(birth_rate) - mll_all - log(alpha));
        delta_j(isinf(delta_j)) = realmax;
        
        delta_all = sum(delta_j);
        
        s = exprnd(1/(birth_rate + delta_all));
        t_fa = t_fa + s;
        
        if(binornd(1,birth_rate/(birth_rate + delta_all)) == 1)
            if theta_b(j).p < p_max
                
%                 disp('birth')
                prior.x0 = zeros(1,1);
                prior.Q0 = eye(1);
                theta_tmp = sample_prior(prior, T, 1, false);
                
                theta_b(j).b = [theta_b(j).b; theta_tmp.b(2)];
                theta_b(j).A = diag([diag(theta_b(j).A);theta_tmp.A(2,2)]);
                theta_b(j).Q = diag([diag(theta_b(j).Q);theta_tmp.Q(2,2)]);
                theta_b(j).X = [theta_b(j).X; theta_tmp.X];
                
                C_b(obsIdx,theta_b(j).p+1) = normrnd(0, ones(N_tmp, 1));
                theta_b(j).p = theta_b(j).p + 1;
            end
        else
%             disp('death')
            delt_j_tmp = delta_j/max(delta_j);
            del_col = mnrnd(1,delt_j_tmp/sum(delt_j_tmp));
            
            del_col_expand = [0 del_col];
            theta_b(j).b = theta_b(j).b(~del_col_expand);
            theta_b(j).A = theta_b(j).A(~del_col_expand, ~del_col_expand);
            theta_b(j).Q = theta_b(j).Q(~del_col_expand, ~del_col_expand);
            theta_b(j).X = theta_b(j).X(~del_col,:);
            
            C_tmp = 0*C_b(obsIdx, :);
            C_tmp(:,1:(theta_b(j).p - 1)) = C_b(obsIdx, ~del_col);
            C_b(obsIdx, :) =  C_tmp;
            theta_b(j).p = theta_b(j).p - 1;
        end
        
    end
    
end




end