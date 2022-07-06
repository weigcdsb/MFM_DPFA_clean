function [theta_b, C_b, delt_b] = constraintProj(theta_a, C_a, delt_a,...
    lab_unique, lab_all)

% to debug
% theta_a = THETA{g};
% C_a = C_fit(:,:,g);
% delt_a = delt_fit(:,g);
% nClus = nClus;
% lab_tmp = Lab;

theta_b = theta_a;
C_b = C_a;
delt_b = delt_a;

for j = lab_unique(:)'
    obsIdx = find(lab_all == j);
    N_tmp = length(obsIdx);
    
    % for mu
    muBar = mean(theta_b(j).mu);
    theta_b(j).mu = theta_b(j).mu - muBar;
    theta_b(j).b(1) = (theta_b(j).A(1,1) - 1)*muBar + theta_b(j).b(1);
    delt_b(obsIdx) = delt_b(obsIdx) + muBar*ones(N_tmp,1);
    
    % for X
    XBar = mean(theta_b(j).X, 2);
    theta_b(j).X = theta_b(j).X - XBar;
    theta_b(j).b(2:end) = (theta_b(j).A(2:end,2:end) - eye(theta_b(j).p))*XBar +...
        theta_b(j).b(2:end);
    delt_b(obsIdx) = delt_b(obsIdx) + C_b(obsIdx,1:theta_b(j).p)*XBar;
end


end