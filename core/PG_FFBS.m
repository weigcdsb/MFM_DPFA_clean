function BETA_b = PG_FFBS(BETA_a, Y_tmp,X_tmp,A,b,...
    Sig,m0,V0,R_tmp,delta_tmp)


% for debugging
% BETA_a = [theta_a.mu;theta_a.X];
% Y_tmp = Y_tmp;
% X_tmp = [ones(N_tmp,1) C_tmp];
% A = theta_a.A;
% b = theta_a.b;
% Sig = theta_a.Q;
% m0 = [prior.mu0;prior.x0];
% V0 = blkdiag(prior.Sigmu0, prior.Q0);
% R_tmp = R_tmp;
% delta_tmp = delta_tmp;

BETA_b = BETA_a;

p = size(BETA_a,1);
T = size(BETA_a,2);
N_tmp = size(Y_tmp,1);

%% (1) calculate r_{nt}, \hat{w}_{nt}, \hat{y}_t
ETA_tmp = X_tmp*BETA_a + delta_tmp;

b_tmp = R_tmp + Y_tmp;
c_tmp = ETA_tmp - log(R_tmp);

b_tmp2 = b_tmp(:);
c_tmp2 = c_tmp(:);
w_tmp_raw = nan*ones(N_tmp*T,1);

% PG draw
obs_idx = ~isnan(b_tmp2);
w_tmp_raw(obs_idx) = ...
pgdraw_expand(b_tmp2(obs_idx), c_tmp2(obs_idx));
w_tmp = reshape(w_tmp_raw, N_tmp,[]);

k_tmp = (Y_tmp - R_tmp)/2 + w_tmp.*...
    (log(R_tmp) - repmat(delta_tmp,1,T));
Yhat_tmp = (1./w_tmp).*k_tmp;

%% (2) FFBS
BETA_b = FFBS(A,b,m0,V0,Sig,Yhat_tmp,X_tmp,w_tmp,T,p);

end

                
                
                
                
                
                