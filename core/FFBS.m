function BETA_b = FFBS(A,b,m0,V0,Sig,...
    Yhat_tmp,X_tmp,w_tmp,T,p)


m_tmp = zeros(p,T);
V_tmp = zeros(p,p,T);

%% (1) FF: forward-filtering
for t = 1:T
    if t == 1
        m_tt_1 = A*m0 + b;
        V_tt_1 = A*V0*A' + Sig;
    else
        m_tt_1 = A*m_tmp(:,t-1) + b;
        V_tt_1 = A*V_tmp(:,:,t-1)*A' + Sig;
    end
    
    obs_idx = ~isnan(Yhat_tmp(:,t));
    X_tmp2 = X_tmp(obs_idx,:);
    w_tmp2 = w_tmp(obs_idx, t);
    Yhat_tmp2 = Yhat_tmp(obs_idx, t);
    
    Kt = V_tt_1*X_tmp2'/(X_tmp2*V_tt_1*X_tmp2' + diag(1./w_tmp2));
    m_tmp(:,t) = m_tt_1 + Kt*(Yhat_tmp2 - X_tmp2*m_tt_1);
    V_tmp(:,:,t) = (eye(p) - Kt*X_tmp2)*V_tt_1;
    
    V_tmp(:,:,t) = (V_tmp(:,:,t) + V_tmp(:,:,t)')/2;
end

%% (2) BS: backward-sampling
BETA_b(:,T) = myMvnrnd(m_tmp(:,T),V_tmp(:,:,T));

for t = (T-1):-1:1
    Jt = V_tmp(:,:,t)*A'/(A*V_tmp(:,:,t)*A' + Sig);
    mstar_tmp = m_tmp(:,t) + Jt*(BETA_b(:,t+1) - A*m_tmp(:,t) - b);
    Vstar_tmp = (eye(p) - Jt*A)*V_tmp(:,:,t);
    Vstar_tmp = (Vstar_tmp + Vstar_tmp')/2;
    BETA_b(:,t) = myMvnrnd(mstar_tmp,Vstar_tmp);
end

end