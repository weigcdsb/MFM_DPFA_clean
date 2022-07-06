function theta_b = update_clusParam_normal(theta_a,delta_tmp,...
    C_tmp,Y_tmp,prior)


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
gradHess = @(vecmuX) gradHessX_na(vecmuX, delta_tmp,...
    [ones(N_tmp,1) C_tmp],...
    [prior.mu0;prior.x0],  blkdiag(prior.Sigmu0, prior.Q0),...
    theta_a.Q, theta_a.A, theta_a.b, Y_tmp);

muX = [theta_a.mu;theta_a.X];
warning('off');
[mudXvec,~,hess_tmp,~, eoi] = newtonGH(gradHess,muX(:),1e-6,1000);
warning('on');

if((sum(isnan(mudXvec)) ~= 0) || eoi == 1000)
    try
        disp('use adaptive smoother initial')
        muX_addp = ppasmoo_poissexp_na(Y_tmp,[ones(N_tmp,1) C_tmp],...
            delta_tmp,...
            [prior.mu0;prior.x0],blkdiag(prior.Sigmu0, prior.Q0),...
            theta_a.A,theta_a.b,theta_a.Q);
        muX_tmp = muX_addp;
    catch
        disp('use 0')
        muX_tmp = muX*0;
    end
    warning('off');
    [mudXvec,~,hess_tmp,~] = newtonGH(gradHess,muX_tmp(:),1e-4,1000);
    warning('on');
end


if sum(isnan(mudXvec)) ~= 0
    muXsamp_trans = muX_addp;
else
    RmuX = chol(-hess_tmp,'lower'); % sparse
    zmuX = randn(length(mudXvec), 1) + RmuX'*mudXvec;
    muXsamp = RmuX'\zmuX;
    muXsamp_trans = reshape(muXsamp,[], T);
end

theta_b.mu = muXsamp_trans(1,:);
theta_b.X = muXsamp_trans(2:end,:);


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