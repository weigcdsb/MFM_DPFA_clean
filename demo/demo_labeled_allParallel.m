usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));

%% generate data
close all
rng(2)
n = 5;
nClus = 2;
N = n*nClus;
p = 2;
T = 1000;

Lab = repelem(1:nClus, n);

delta = randn(n*nClus,1)*0.5;
C_trans = zeros(n*nClus, p*nClus);
for k = 1:nClus
    NtmpIdx = (n*(k-1)+1):(n*k);
    delta(NtmpIdx) = delta(NtmpIdx) - mean(delta(NtmpIdx));
    pIdx = id2id(k, p);
    C_trans(NtmpIdx, pIdx) = randn(n, p);
end


X = ones(p*nClus, T)*Inf;
x0 = zeros(p*nClus,1);
Q0 = eye(nClus*p)*1e-2;
X(:,1) = mvnrnd(x0, Q0)';

% Generate X offline (A unspecified)
for i=1:(nClus*p)
    while(sum(X(i,:) > 1.5) > 1)
        k = ceil(rand()*20)+10;
        X(i,:) = interp1(linspace(0,1,k),randn(k,1)*0.2,linspace(0,1,T),'spline');
    end
end

% generate mu
mu = ones(nClus, T)*Inf;
for i = 1:nClus
    while(sum(mu(i,:) > log(15)) > 1 || sum(mu(i,:) < log(2)) == T)
        k = ceil(rand()*25)+10;
        mu(i,:) = interp1(linspace(0,1,k),randn(k,1)*0.7,linspace(0,1,T),'spline');
    end
end

mu = mu - mean(mu,2);

% let's generate lambda
logLam = zeros(n*nClus, T);
logLam(:,1) = delta + C_trans*X(:,1) + repelem(mu(:,1), n);


for t=2:T
    logLam(:, t) = delta + C_trans*X(:,t) + repelem(mu(:,t), n);
end
%
Y = poissrnd(exp(logLam));

% extract data out?
for kk=1:N
    Y(kk,randsample(T,round(T/2))) = nan;
end

%% pre-MCMC
rng(1)
ng = 1000;
burnIn = round(ng/20);

p_max = 20;

% priors...
prior.mu0 = 0;
prior.Sigmu0 = 1;

prior.delt0 = 0;
prior.Sigdelt0 = 1;

prior.BA0 =[0 1]';
prior.Lamb0 = eye(2);
prior.nu0 = 1;
prior.sig20 = 1e-2;

% pre-allocation
delt_fit = zeros(N,ng);
C_fit = zeros(N,p_max,ng);

delt_fit(:,1) = normrnd(prior.delt0, prior.Sigdelt0, N,1);
p_trace = ones(nClus, ng);

for k = 1:nClus
    idx_tmp = find(Lab == k);
    N_tmp = length(idx_tmp);
    p_tmp = p_trace(k,1);
    
    prior.x0 = zeros(p_tmp,1);
    prior.Q0 = eye(p_tmp);
    prior.muC0 = zeros(p_tmp,1);
    prior.SigC0 = eye(p_tmp);
    
    C_fit(idx_tmp,1:p_tmp, 1) = mvnrnd(prior.muC0, prior.SigC0, N_tmp);
    theta_tmp = sample_prior(prior, T, p_tmp, false);
    theta_tmp.p = p_tmp;
    
    THETA{1}(k) = theta_tmp;
end

% llhd_trace
fitMFRTrace = zeros(N, T, ng);
llhd_trace = zeros(ng,1);

% tuning parameter: NB dispersion
R_all = 10*ones(N,T);

% acceptance rate for PG-FFBS-MH
ACC_trace = zeros(nClus, ng);
RACC = zeros(nClus, 1);

% BDMCMC
birth_rate = .5;
birth_time = 1;
alpha = 2;

% AGS
bdstart = 50;
b0 = 0.1;
b1 = 1e-3;

%% MCMC
% poolobj = parpool('IdleTimeout',120);

tic
for g= 2:ng
    
    % (1) update cluster parameters (PG+FFBS+MH)
    
    
    THETA_tmp_in = THETA{g-1};
    THETA_tmp_out = THETA{g-1};
    ACC_tmp = ACC_trace(:,g);
    
    parfor j = 1:nClus
        obsIdx = find(Lab == j);
        
        prior_tmp = prior;
        prior_tmp.x0 = zeros(THETA_tmp_in(j).p,1);
        prior_tmp.Q0 = eye(THETA_tmp_in(j).p);
        prior_tmp.muC0 = zeros(THETA_tmp_in(j).p,1);
        prior_tmp.SigC0 = eye(THETA_tmp_in(j).p);
        
        [THETA_tmp_out(j), acc] = update_clusParam_PG(THETA_tmp_in(j),...
            delt_fit(obsIdx,g-1), C_fit(obsIdx,1:THETA_tmp_in(j).p,g-1),...
            Y(obsIdx,:),R_all(obsIdx,:),prior_tmp);
        ACC_tmp(j) = acc;
    end
    
    THETA{g} = THETA_tmp_out;
    ACC_trace(:,g) = ACC_tmp;
    
    for j = 1:nClus
        if(g > burnIn)
            RACC(j) = sum(ACC_trace(j, (burnIn+1):g))/(g - burnIn);
        end
    end
    
    % (2) update non-cluster parameters: intercept + loading (HMC)
    [delt_fit(:,g), C_fit(:,:,g)] =...
        update_deltC_HMC(delt_fit(:,g-1), C_fit(:,:,g-1),...
        THETA{g},Lab,Y,prior);
    
    % (3) projection
    for j = 1:nClus
        obsIdx = find(Lab == j);
        N_tmp = length(obsIdx);
        
        % for mu
        muBar = mean(THETA{g}(j).mu);
        THETA{g}(j).mu = THETA{g}(j).mu - muBar;
        THETA{g}(j).b(1) = (THETA{g}(j).A(1,1) - 1)*muBar + THETA{g}(j).b(1);
        delt_fit(obsIdx,g) = delt_fit(obsIdx,g) + muBar*ones(N_tmp,1);
        
        % for X
        XBar = mean(THETA{g}(j).X, 2);
        THETA{g}(j).X = THETA{g}(j).X - XBar;
        THETA{g}(j).b(2:end) = (THETA{g}(j).A(2:end,2:end) - eye(THETA{g}(j).p))*XBar +...
            THETA{g}(j).b(2:end);
        delt_fit(obsIdx,g) = delt_fit(obsIdx,g) + C_fit(obsIdx,1:THETA{g}(j).p,g)*XBar;
    end
    
    fitMFR = zeros(N, T);
    for k  = 1:N
        fitMFR(k,:) = exp([1 delt_fit(k,g) C_fit(k,1:THETA{g}(Lab(k)).p,g)]*...
            [THETA{g}(Lab(k)).mu ;ones(1,T) ;THETA{g}(Lab(k)).X]);
    end
    fitMFRTrace(:,:,g) = fitMFR;
    llhd_trace(g) = nansum(log(poisspdf(Y,fitMFR)), 'all')/...
        nansum(Y, 'all');
    
    
    % birth-death
    if g > bdstart && binornd(1, exp(-b0 - b1*g)) == 1
        
        THETA_tmp = THETA{g};
        C_tmp_cell = cell(nClus,1);
        
        for j = 1:nClus
            obsIdx = find(Lab == j);
            C_tmp_cell{j} = C_fit(obsIdx,:,g);
        end
        
        
        parfor j = 1:nClus
            
            t_fa = 0;
            obsIdx = find(Lab == j);
            N_tmp = length(obsIdx);
            
            while(t_fa < birth_time)
                
                mll_del = zeros(THETA_tmp(j).p,1);
                if THETA_tmp(j).p > 1
                    for kk = 1:THETA_tmp(j).p
                        colSelect = setdiff(1:THETA_tmp(j).p, kk);
                        mll_del(kk) = poiLogMarg(reshape(Y(obsIdx,:)', [], 1),...
                            repmat(THETA_tmp(j).X(colSelect,:)', N_tmp,1),...
                            kron(delt_fit(obsIdx,g), ones(T,1)) +...
                            repmat(THETA_tmp(j).mu', N_tmp,1));
                    end
                    
                    mll_all = poiLogMarg(reshape(Y(obsIdx,:)', [], 1),...
                        repmat(THETA_tmp(j).X', N_tmp,1),...
                        kron(delt_fit(obsIdx,g), ones(T,1)) +...
                        repmat(THETA_tmp(j).mu', N_tmp,1));
                else
                    for kk = 1:THETA_tmp(j).p
                        colSelect = setdiff(1:THETA_tmp(j).p, kk);
                        
                        eta_tmp =...
                            [ones(N_tmp,1) delt_fit(obsIdx,g) C_tmp_cell{j}(:,colSelect)]*...
                            [THETA_tmp(j).mu ;ones(1,T) ;THETA_tmp(j).X(colSelect,:)];
                        mll_del(kk) = nansum(log(poisspdf(Y(obsIdx,:),exp(eta_tmp))), 'all');
                    end
                    
                    eta_tmp =...
                        [ones(N_tmp,1) delt_fit(obsIdx,g) C_tmp_cell{j}(:,1:THETA_tmp(j).p)]*...
                        [THETA_tmp(j).mu ;ones(1,T) ;THETA_tmp(j).X];
                    
                    mll_all = nansum(log(poisspdf(Y(obsIdx,:),exp(eta_tmp))), 'all');
                end
                
                delta_j = exp(mll_del + log(birth_rate) - mll_all - log(alpha));
                delta_j(isinf(delta_j)) = realmax;
                
                delta_all = sum(delta_j);
                s = exprnd(1/(birth_rate + delta_all));
                t_fa = t_fa + s;
                
                if(binornd(1,birth_rate/(birth_rate + delta_all)) == 1)
                    if THETA_tmp(j).p < p_max
                        % give a birth
                        disp('birth')
                        prior_tmp = prior;
                        
                        prior_tmp.x0 = zeros(1,1);
                        prior_tmp.Q0 = eye(1);
                        theta_tmp = sample_prior(prior_tmp, T, 1, false);
                        
                        THETA_tmp(j).b = [THETA_tmp(j).b; theta_tmp.b(2)];
                        THETA_tmp(j).A = diag([diag(THETA_tmp(j).A);theta_tmp.A(2,2)]);
                        THETA_tmp(j).Q = diag([diag(THETA_tmp(j).Q);theta_tmp.Q(2,2)]);
                        THETA_tmp(j).X = [THETA_tmp(j).X; theta_tmp.X];
                        
                        C_tmp_cell{j}(:,THETA_tmp(j).p+1) = normrnd(0,1, N_tmp,1);
                        THETA_tmp(j).p = THETA_tmp(j).p + 1;
                    end
                else
                    
                    % give a death
                    disp('death')
                    del_col = mnrnd(1,delta_j/delta_all);
                    
                    del_col_expand = [0 del_col];
                    THETA_tmp(j).b = THETA_tmp(j).b(~del_col_expand);
                    THETA_tmp(j).A = THETA_tmp(j).A(~del_col_expand, ~del_col_expand);
                    THETA_tmp(j).Q = THETA_tmp(j).Q(~del_col_expand, ~del_col_expand);
                    THETA_tmp(j).X = THETA_tmp(j).X(~del_col,:);
                    
                    C_tmp = 0*C_tmp_cell{j};
                    C_tmp(:,1:(THETA_tmp(j).p - 1)) = C_tmp_cell{j}(:,~del_col);
                    
                    C_tmp_cell{j} =  C_tmp;
                    THETA_tmp(j).p = THETA_tmp(j).p - 1;
                end
                
            end
        end
        
        THETA{g} = THETA_tmp;
        for j = 1:nClus
            obsIdx = find(Lab == j);
            C_fit(obsIdx,:,g) = C_tmp_cell{j};
        end
    end
    
    p_trace(1,g) = THETA{g}(1).p;
    p_trace(2,g) = THETA{g}(2).p;
    
    if g > bdstart
        figure(1)
        subplot(1,2,1)
        histogram(p_trace(1,bdstart:g))
        xline(2, 'r--', 'Linewidth', 2)
        title("p^{(g)} = "+p_trace(1,g))
        subplot(1,2,2)
        histogram(p_trace(2,bdstart:g))
        xline(2, 'r--', 'Linewidth', 2)
        title("p^{(g)} = "+p_trace(2,g))
    end
    
    
    figure(2)
    hold on
    plot(delt_fit(:,g), 'k')
    plot(delta, 'r')
    hold off
    
    figure(3)
    subplot(2,2,1)
    plot(THETA{g}(1).mu)
    ylabel("racc = " + RACC(1))
    subplot(2,2,2)
    plot(mu(1,:))
    
    subplot(2,2,3)
    plot(THETA{g}(2).mu)
    ylabel("racc = " + RACC(2))
    subplot(2,2,4)
    plot(mu(2,:))
    
    disp("iter " + g +": "+ toc)
end
toc








