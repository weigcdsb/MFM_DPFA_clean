function [Z_b, numClus_b, t_b, actList_b,...
    THETA_b,delt_b,C_b] =...
    splitMerge(Z_a, numClus_a, t_a, actList_a,...
    THETA_a,delt_a,C_a,Y_tmp,R_all, prior,...
    a, b, log_v,splitL, nSplit, nMerge)


% to debug
% Z_a = Z_fit(:,g);
% numClus_a = numClus_fit(:,g);
% t_a = t_fit(g);
% actList_a = actList(:,g);
% c_next_a = ordered_next(actList(:,g));
% THETA_a = THETA{g};
% delt_a = delt_fit(:,g);
% C_a = C_fit(:,:,g);
% splitL = 1;
% Y_tmp = Y;
% R_all = R_all;
% nSplit = 5;
% nMerge = 5;


Z_b = Z_a;
numClus_b = numClus_a;
t_b = t_a;
actList_b = actList_a;
THETA_b = THETA_a;
delt_b = delt_a;
C_b = C_a;

N = size(Y_tmp,1);
T = size(Y_tmp,2);
lAbsGam = @(x) log(abs(gamma(x)));

% (a) randomly choose a pair of indices
if ~isnan(splitL)
    if splitL
        [~,maxIdx] = max(numClus_b);
        rdIdx = randsample(find(Z_b == maxIdx),2);
    else
        sVal2 = mink(numClus_b(numClus_b ~= 0), 2);
        
        if length(sVal2) == 1
            rdIdx = randsample(N,2);
        else
            set1 = find(numClus_b == sVal2(1));
            minIdx1 = set1(randsample(length(set1),1));
            
            set2 = setdiff(find(numClus_b == sVal2(2)), minIdx1);
            minIdx2 = set2(randsample(length(set2),1));
            
            
            if length(find(Z_b == minIdx1)) > 1
                rdIdx = randsample(find(Z_b == minIdx1),1);
            else
                rdIdx = find(Z_b == minIdx1);
            end
            
            if length(find(Z_b == minIdx2)) > 1
                rdIdx = [rdIdx; randsample(find(Z_b == minIdx2),1)];
            else
                rdIdx = [rdIdx; find(Z_b == minIdx2)];
            end
            
        end
    end
else
    rdIdx = randsample(N,2);
end


% [~,maxIdx] = max([THETA_b(Z_b(rdIdx(1))).p THETA_b(Z_b(rdIdx(2))).p]);

ism = rdIdx(1);% rdIdx(maxIdx);
jsm = rdIdx(2);% rdIdx(setdiff(1:2, maxIdx));

ci0 = Z_b(ism);
cj0 = Z_b(jsm);
ti0 = THETA_b(ci0);
tj0 = THETA_b(cj0);


% (b) set S(1),...,S(ns) to the indices of
% the points in clusters ci0 and cj0
S = zeros(N,1);
ns = 0;
for k = 1:N
    if(Z_b(k) == ci0 || Z_b(k) == cj0)
        ns = ns + 1;
        S(ns) = k;
    end
end

% (c) find available cluster IDs
% for merge and split parameters
k = 1;
while(actList_b(k) == k); k = k+1;end;cm = k;
while(actList_b(k) == k+1); k = k+1;end;ci = k+1;
while(actList_b(k) == k+2); k = k+1;end;cj = k+2;

% (d) randomly choose the merge launch state
psm_min = min(THETA_a(ci0).p, THETA_a(cj0).p);
prior.x0 = zeros(psm_min,1);
prior.Q0 = eye(psm_min);
prior.muC0 = zeros(psm_min,1);
prior.SigC0 = eye(psm_min);

tm = sample_prior(prior, T, psm_min, false); % ti0
tm.p = psm_min;
mIdx = S(1:ns);

delt_m = delt_b;
C_m = C_b;
C_m(mIdx, (tm.p+1):end) = 0;


% useNormal = false;
for mm = 1:nMerge
    
    % update factor
    if mm<2
        useNormal = true;
    else
        useNormal = false;
    end
    
    try
        [tm, ~] = update_factor_sm(tm,tm,delt_m(mIdx),...
            C_m(mIdx,1:tm.p),Y_tmp(mIdx,:),...
            prior,true,false, R_all(mIdx,:),useNormal);
        
    catch
        [tm, ~] = update_factor_sm(tm,tm,delt_m(mIdx),...
            C_m(mIdx,1:tm.p),Y_tmp(mIdx,:),...
            prior,true,false, R_all(mIdx,:),false);
    end
    
    
    
    % update loading
    [delt_m(mIdx), C_m(mIdx,:)] =...
        update_deltC_HMC(delt_m(mIdx), C_m(mIdx,:),...
        tm,ones(size(mIdx)),Y_tmp(mIdx,:),prior);
    
end


% (e) randomly choose the split lauch state
prior.x0 = zeros(ti0.p,1);
prior.Q0 = eye(ti0.p);
prior.muC0 = zeros(ti0.p,1);
prior.SigC0 = eye(ti0.p);
ti = sample_prior(prior, T, ti0.p, false); % ti0
ti.p = ti0.p;

prior.x0 = zeros(tj0.p,1);
prior.Q0 = eye(tj0.p);
prior.muC0 = zeros(tj0.p,1);
prior.SigC0 = eye(tj0.p);
tj = sample_prior(prior, T, tj0.p, false); % tj0
tj.p = tj0.p;

zs = ones(N,1);

zs(ism) = ci;
zs(jsm) = cj;

kOut = setdiff(S(1:ns), [ism,jsm]);
splitIdx = binornd(1,.5,length(kOut), 1);
siIdx = kOut(splitIdx == 1);
sjIdx = kOut(splitIdx == 0);
zs(siIdx) = ci;
zs(sjIdx) = cj;
iIdx = [ism;siIdx];
jIdx = [jsm;sjIdx];
ni = length(iIdx);
nj = length(jIdx);

% make several moves (restricted Gibbs)
delt_s = delt_b;
C_s = C_b;

C_s(iIdx, (ti.p+1):end) = 0;
C_s(jIdx, (tj.p+1):end) = 0;

% useNormal = false;
for ss = 1:nSplit
    
    if ss<2
        useNormal = true;
    else
        useNormal = false;
    end
    
    % update factors...
    try
        [zs,ti,tj,ci,cj,~, ni, nj, iIdx, jIdx] =...
            restricted_gibbs_sm(zs,zs,ti,ti,tj,tj,ci,ci,cj,cj,...
            iIdx,jIdx, ni,nj,ism,jsm,S,ns,Y_tmp,b,prior, true, delt_s, C_s, R_all,useNormal);
    catch
        [zs,ti,tj,ci,cj,~, ni, nj, iIdx, jIdx] =...
            restricted_gibbs_sm(zs,zs,ti,ti,tj,tj,ci,ci,cj,cj,...
            iIdx,jIdx, ni,nj,ism,jsm,S,ns,Y_tmp,b,prior, true, delt_s, C_s, R_all,false);
    end
    
    
    C_s(iIdx, (ti.p+1):end) = 0;
    C_s(jIdx, (tj.p+1):end) = 0;
    
    % update loadings...
    [delt_s(iIdx), C_s(iIdx,:)] =...
        update_deltC_HMC(delt_s(iIdx), C_s(iIdx,:), ti,ones(size(iIdx)),...
        Y_tmp(iIdx,:),prior);
    [delt_s(jIdx), C_s(jIdx,:)] =...
        update_deltC_HMC(delt_s(jIdx), C_s(jIdx,:), tj,ones(size(jIdx)),...
        Y_tmp(jIdx,:),prior);
end


% (f) make proposal
if ci0 == cj0 % propose a split
    % make one final sweep and compute its density
    % update factors...
    
    [zs,ti,tj,ci,cj,log_prop_ab, ni, nj, iIdx, jIdx] =...
        restricted_gibbs_sm(zs,zs,ti,ti,tj,tj,ci,ci,cj,cj,...
        iIdx,jIdx, ni,nj,ism,jsm,S,ns,Y_tmp,b,prior, true, delt_s, C_s, R_all,false);
    C_s(iIdx, (ti.p+1):end) = 0;
    C_s(jIdx, (tj.p+1):end) = 0;
    
    % update loadings...
    [delt_s(iIdx), C_s(iIdx,:)] =...
        update_deltC_HMC(delt_s(iIdx), C_s(iIdx,:), ti,ones(size(iIdx)),...
        Y_tmp(iIdx,:),prior);
    [delt_s(jIdx), C_s(jIdx,:)] =...
        update_deltC_HMC(delt_s(jIdx), C_s(jIdx,:), tj,ones(size(jIdx)),...
        Y_tmp(jIdx,:),prior);
    
    % compute density of Lmerge to original state
    [ti0, log_prop_ba] = update_factor_sm(tm,ti0,delt_a(mIdx),...
        C_a(mIdx,1:ti0.p),Y_tmp(mIdx,:),...
        prior,false,true, R_all(mIdx,:),false);
    
    % compute acceptance probability
    log_prior_b = log_v(t_b+1) + lAbsGam(ni+b) +...
        lAbsGam(nj+b) - 2*lAbsGam(a) +...
        log_prior_sm(ti, prior, delt_s(iIdx,:), C_s(iIdx,:));
    log_prior_a = log_v(t_b) + lAbsGam(ns+b) - lAbsGam(a) +...
        log_prior_sm(ti0, prior, delt_s(Z_b == ci0,:), C_s(Z_b == ci0,:));
    
    llhd_ratio = 0;
    for ks = 1:ns
        k = S(ks);
        if zs(k) == ci
            lamTmp = exp([1 delt_s(k) C_s(k,1:ti.p)]*...
                [ti.mu ;ones(1,T) ;ti.X]);
        else
            lamTmp = exp([1 delt_s(k) C_s(k,1:tj.p)]*...
                [tj.mu ;ones(1,T) ;tj.X]);
        end
        
        lamTmp0 = exp([1 delt_a(k) C_a(k,1:ti0.p)]*...
            [ti0.mu ;ones(1,T) ;ti0.X]);
        llhd_ratio = llhd_ratio +...
            nansum(log(poisspdf(Y_tmp(k,:), lamTmp))) -...
            nansum(log(poisspdf(Y_tmp(k,:), lamTmp0)));
    end
    p_accept = min(1, exp(log_prop_ba-log_prop_ab +...
        log_prior_b-log_prior_a + llhd_ratio));
    
    if rand < p_accept % accept split
        disp('accept split')
        for ks = 1:ns
            Z_b(S(ks)) = zs(S(ks));
        end
        actList_b = ordered_remove(ci0, actList_b, t_b);
        actList_b = ordered_insert(ci, actList_b, t_b-1);
        actList_b = ordered_insert(cj, actList_b, t_b);
        
        numClus_b(ci0) = 0;
        numClus_b(ci) = ni;
        numClus_b(cj) = nj;
        t_b = t_b + 1;
        
        THETA_b(ci) = ti;
        THETA_b(cj) = tj;
        delt_b = delt_s;
        C_b = C_s;
        
        idxAll = [ci cj];
        for j = 1:2
            c = idxAll(j);
            obsIdx = find(Z_b == c);
            N_tmp = length(obsIdx);
            
            % for mu
            muBar = mean(THETA_b(c).mu);
            THETA_b(c).mu = THETA_b(c).mu - muBar;
            THETA_b(c).b(1) = (THETA_b(c).A(1,1) - 1)*muBar + THETA_b(c).b(1);
            delt_b(obsIdx) = delt_b(obsIdx) + muBar*ones(N_tmp,1);
            
            % for X
            XBar = mean(THETA_b(c).X, 2);
            THETA_b(c).X = THETA_b(c).X - XBar;
            THETA_b(c).b(2:end) = (THETA_b(c).A(2:end,2:end) - eye(psm_min))*XBar +...
                THETA_b(c).b(2:end);
            delt_b(obsIdx) = delt_b(obsIdx) + C_b(obsIdx,1:THETA_b(c).p)*XBar;
        end
    else
        disp('reject split')
    end
else % propose a merge
    % make one final sweep and compute its probability density
    
    % update factor
    [tm, log_prop_ab] = update_factor_sm(tm,tm,delt_m(mIdx),...
        C_m(mIdx,1:tm.p),Y_tmp(mIdx,:),...
        prior,true,true, R_all(mIdx,:),false);
    
    % update loading
    [delt_m(mIdx), C_m(mIdx,:)] =...
        update_deltC_HMC(delt_m(mIdx), C_m(mIdx,:),...
        tm,ones(size(mIdx)),Y_tmp(mIdx,:),prior);
    
    % compute probability density of going from split launch state to original state
    [~,~,~,~,~,log_prop_ba, ni, nj, iIdx, jIdx] =...
        restricted_gibbs_sm(zs,Z_b,ti,ti0,tj,tj0,ci,ci0,cj,cj0,...
        iIdx,jIdx, ni,nj,ism,jsm,S,ns,Y_tmp,b,prior, false, delt_s, C_s, R_all,false);
    
    % compute acceptance probability
    log_prior_b = log_v(t_b-1) + lAbsGam(ns+b)-lAbsGam(a) +...
        log_prior_sm(tm, prior, delt_s(mIdx,:), C_s(mIdx,:));
    log_prior_a = log_v(t_b) + lAbsGam(ni+b)+lAbsGam(nj+b)-2*lAbsGam(a) +...
        log_prior_sm(ti0, prior, delt_s(Z_b == ci0,:), C_s(Z_b == ci0,:)) +...
        log_prior_sm(tj0, prior, delt_s(Z_b == cj0,:), C_s(Z_b == cj0,:));
    llhd_ratio = 0;
    for ks = 1:ns
        k = S(ks);
        
        if Z_b(k) == ci0
            lamTmp0 = exp([1 delt_a(k) C_a(k,1:ti0.p)]*...
                [ti0.mu ;ones(1,T) ;ti0.X]);
        else
            lamTmp0 = exp([1 delt_a(k) C_a(k,1:tj0.p)]*...
                [tj0.mu ;ones(1,T) ;tj0.X]);
        end
        lamTmp = exp([1 delt_m(k) C_m(k,1:tm.p)]*...
            [tm.mu ;ones(1,T) ;tm.X]);
        
        llhd_ratio = llhd_ratio +...
            nansum(log(poisspdf(Y_tmp(k,:), lamTmp))) -...
            nansum(log(poisspdf(Y_tmp(k,:), lamTmp0)));
    end
    p_accept = min(1.0, exp(log_prop_ba-log_prop_ab +...
        log_prior_b-log_prior_a + llhd_ratio));
    
    if rand < p_accept % accept merge
        disp('accept merge')
        for ks = 1:ns
            Z_b(S(ks)) = cm;
        end
        actList_b = ordered_remove(ci0, actList_b, t_b);
        actList_b = ordered_remove(cj0, actList_b, t_b-1);
        actList_b = ordered_insert(cm, actList_b, t_b-2);
        
        numClus_b(cm) = ns;
        numClus_b(ci0) = 0;
        numClus_b(cj0) = 0;
        t_b = t_b - 1;
        
        THETA_b(cm) = tm;
        delt_b = delt_m;
        C_b = C_m;
        
        % projection
        obsIdx = find(Z_b == cm);
        N_tmp = length(obsIdx);
        
        % for mu
        muBar = mean(THETA_b(cm).mu);
        THETA_b(cm).mu = THETA_b(cm).mu - muBar;
        THETA_b(cm).b(1) = (THETA_b(cm).A(1,1) - 1)*muBar + THETA_b(cm).b(1);
        delt_b(obsIdx) = delt_b(obsIdx) + muBar*ones(N_tmp,1);
        
        % for X
        XBar = mean(THETA_b(cm).X, 2);
        THETA_b(cm).X = THETA_b(cm).X - XBar;
        THETA_b(cm).b(2:end) = (THETA_b(cm).A(2:end,2:end) - eye(psm_min))*XBar +...
            THETA_b(cm).b(2:end);
        delt_b(obsIdx) = delt_b(obsIdx) + C_b(obsIdx,1:THETA_b(cm).p)*XBar;
        
    else
        disp('reject merge')
    end
end



end