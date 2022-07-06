function [theta_b, log_pdf] = update_factor_sm(theta_a,theta_b,delta_tmp, C_tmp,Y_tmp,...
    prior,active,density, R_tmp, useNormal)

% to debug
% theta_a = tm;
% theta_b = tm;
% delta_tmp = delt_m(mIdx);
% C_tmp = C_m(mIdx,1:theta_a.p);
% Y_tmp = Y_tmp(mIdx,:);
% active =true;
% density = true;
% R_tmp = R_all(mIdx,:);

log_pdf = 0;

N_tmp = size(Y_tmp, 1);
T = size(Y_tmp,2);
p = size(C_tmp, 2);

prior.x0 = zeros(p,1);
prior.Q0 = eye(p);
prior.muC0 = zeros(p,1);
prior.SigC0 = eye(p);


%% (1) update X & mu

if active
    if useNormal
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
    else
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
        else
            theta_b.mu = theta_a.mu;
            theta_b.X = theta_a.X;
        end
    end
end


% update linear dynamics
muX = [theta_b.mu;theta_b.X];

for k = 1:(p+1)
    
    Y_BA = muX(k,2:T)';
    X_BA = [ones(T-1,1) muX(k,1:(T-1))'];
    
    Lamn = X_BA'*X_BA + prior.Lamb0;
    BAn = Lamn\(prior.Lamb0*prior.BA0 + X_BA'*Y_BA);
    
    an = (prior.nu0 + T-1)/2;
    bn = (prior.nu0*prior.sig20)/2 +...
        (Y_BA'*Y_BA + ...
        prior.BA0'*prior.Lamb0*prior.BA0 -...
        BAn'*Lamn*BAn)/2;
    
    if active;theta_b.Q(k,k) = 1/gamrnd(an,1/bn);end
    theta_b.Q(k,k) = min(theta_b.Q(k,k), 1e-2);
    if density
        log_pdf = log_pdf + log(gampdf(1/theta_b.Q(k,k), an, 1/bn));
    end
    
    
    % (6) update b_fit & A_fit
    if active
        BAsamp = myMvnrnd(BAn, theta_b.Q(k,k)*inv(Lamn));
        theta_b.b(k) = BAsamp(1);
        theta_b.A(k,k) = BAsamp(2);
    end
    
    if density
        baTmp = [theta_b.b(k) theta_b.A(k,k)]';
        log_pdf = log_pdf + mvnlpdf(baTmp, BAn(:), theta_b.Q(k,k)*inv(Lamn));
    end
    
end



end