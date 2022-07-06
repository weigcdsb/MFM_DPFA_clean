function [theta_b, acc_b] = update_clusParam_all(theta_a, C_a, delt_a, lab_unique,...
    lab_all, Y_tmp, prior, R_all, useNormal)


% to debug
% theta_a = THETA{g-1};
% C_a = C_fit(:,:,g-1);
% delt_a = delt_fit(:,g-1);
% lab_unique = unique(Lab);
% lab_all = Lab;
% Y_tmp = Y;
% prior = prior;
% R_all = R_all;
% useNormal = false;

theta_b = theta_a;
acc_b = zeros(max(lab_unique), 1);

for j = lab_unique(:)'
    obsIdx = find(lab_all == j);
    
    prior_tmp = prior;
    prior_tmp.x0 = zeros(theta_a(j).p,1);
    prior_tmp.Q0 = eye(theta_a(j).p);
    prior_tmp.muC0 = zeros(theta_a(j).p,1);
    prior_tmp.SigC0 = eye(theta_a(j).p);
    
    if useNormal
        theta_b(j) = update_clusParam_normal(theta_a(j),...
            delt_a(obsIdx), C_a(obsIdx,1:theta_a(j).p),...
            Y_tmp(obsIdx,:),prior_tmp);
    else
        [theta_b(j), acc_b(j)] = update_clusParam_PG(theta_a(j),...
            delt_a(obsIdx), C_a(obsIdx,1:theta_a(j).p),...
            Y_tmp(obsIdx,:),R_all(obsIdx,:),prior_tmp);
    end
    
end


end