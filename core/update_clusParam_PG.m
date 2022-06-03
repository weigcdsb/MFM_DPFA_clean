function [theta_b, acc] = update_clusParam_PG(theta_a,delta_tmp,...
    C_tmp,Y_tmp,R_tmp,prior)


% for debugging...
% theta_a = THETA{g-1}(j);
% delta_tmp = delt_fit(obsIdx,g-1);
% C_tmp = C_fit(obsIdx,1:THETA{g-1}(j).p,g-1);
% Y_tmp = Y(obsIdx,:);
% R_tmp = R_all(obsIdx,:);
% prior = prior;


theta_b = theta_a;

N_tmp = size(Y_tmp, 1);
T = size(Y_tmp,2);
p = size(C_tmp, 2);

%% (1) sample proposal
muX_new = PG_FFBS([theta_a.mu;theta_a.X],Y_tmp,...
    [ones(N_tmp,1) C_tmp],theta_a.A,...
    theta_a.b,theta_a.Q,[prior.mu0;prior.x0],...
    blkdiag(prior.Sigmu0, prior.Q0),R_tmp,delta_tmp);


%% (2) MH step
eta_ori = [ones(N_tmp,1) delta_tmp C_tmp]*[theta_a.mu ;ones(1,T); theta_a.X];
eta_new = [ones(N_tmp,1) delta_tmp C_tmp]*...
    [muX_new(1,:) ;ones(1,T); muX_new(2:end,:)];

lam_ori = exp(eta_ori);
lam_new = exp(eta_new);
p_ori = 1./(1 + exp(-eta_ori + log(R_tmp)));
p_new = 1./(1 + exp(-eta_new + log(R_tmp)));

lhr = nansum(log(poisspdf(Y_tmp,lam_new)), 'all')-...
        nansum(log(poisspdf(Y_tmp,lam_ori)), 'all')+...
        nansum(log(nbinpdf(Y_tmp,R_tmp,p_ori)), 'all')-...
        nansum(log(nbinpdf(Y_tmp,R_tmp,p_new)), 'all');

if(log(rand) < lhr)
    theta_b.mu = muX_new(1,:);
    theta_b.X = muX_new(2:end,:);
    acc = 1;
else
    theta_b.mu = theta_a.mu;
    theta_b.X = theta_a.X;
    acc = 0;
end

%% (3) update linear dynamics
muX = [theta_b.mu;theta_b.X];

for k = 1:(p+1)
    try
        Y_BA = muX(k,2:T)';
        X_BA = [ones(T-1,1) muX(k,1:(T-1))'];
        
        Lamn = X_BA'*X_BA + prior.Lamb0;
        BAn = Lamn\(prior.Lamb0*prior.BA0 + X_BA'*Y_BA);
        
        an = (prior.nu0 + T-1)/2;
        bn = (prior.nu0*prior.sig20)/2 +...
            (Y_BA'*Y_BA + ...
            prior.BA0'*prior.Lamb0*prior.BA0 -...
            BAn'*Lamn*BAn)/2;
        theta_b.Q(k,k) = 1/gamrnd(an,1/bn);
    catch
        theta_b.Q(k,k) = theta_a.Q(k,k);
    end
    
    theta_b.Q(k,k) = min(theta_b.Q(k,k), 1e-1);
    
    % (6) update b_fit & A_fit
    try
        BAsamp = myMvnrnd(BAn, theta_b.Q(k,k)*inv(Lamn));
        theta_b.b(k) = BAsamp(1);
        theta_b.A(k,k) = BAsamp(2);
    catch
        theta_b.b(k) = theta_a.b(k);
        theta_b.A(k,k) = theta_a.A(k,k);
    end
end



end