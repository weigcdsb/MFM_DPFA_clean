usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));

%% generate data
close all
rng(2)
n = 5;
nClus = 1;
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
logLam(:,1) = delta + C_trans*X(:,1) + repelem(mu(:,1), n)';


for t=2:T
    logLam(:, t) = delta + C_trans*X(:,t) + repelem(mu(:,t), n)';
end
%
Y = poissrnd(exp(logLam));

% extract data out?
% for kk=1:N
%     Y(kk,randsample(T,round(T/2))) = nan;
% end

%% pre-MCMC
rng(1)
close all
ng = 10000;
burnIn = round(ng/10);

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
birth_rate = .4;
birth_time = 10*(1/birth_rate);
alpha = 5;

% AGS
bdstart = 2;
b0 = 0.1;
b1 = 1e-4;

%% MCMC
if isempty(gcp('nocreate'))
    poolobj = parpool('IdleTimeout',120);
end


% initialize with Laplace approximation
gPre = 10;
for g = 2:gPre
    
    disp("g = " + g)
    
    % (1) update cluster parameters
    [THETA{1}, ~] = update_clusParam_all(THETA{1}, C_fit(:,:,1),...
        delt_fit(:,1), unique(Lab),...
        Lab, Y, prior, R_all, true);
    
    % (2) update non-cluster parameters: intercept + loading (HMC)
    [delt_fit(:,1), C_fit(:,:,1)] =...
        update_deltC_HMC(delt_fit(:,1), C_fit(:,:,1),...
        THETA{1},Lab,Y,prior);
    
    % (3) projection
    [THETA{1},  C_fit(:,:,1), delt_fit(:,1)] = constraintProj(THETA{1},...
         C_fit(:,:,1), delt_fit(:,1),unique(Lab), Lab);
end

tic
for g= 2:ng
    
    % (1) update cluster parameters (PG+FFBS+MH)
    [THETA{g}, ACC_trace(:,g)] = update_clusParam_all(THETA{g-1}, C_fit(:,:,g-1),...
        delt_fit(:,g-1), unique(Lab),...
        Lab, Y, prior, R_all, false);
    
    
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
    [THETA{g},  C_fit(:,:,g), delt_fit(:,g)] = constraintProj(THETA{g},...
         C_fit(:,:,g), delt_fit(:,g),unique(Lab), Lab);
    
    
    fitMFR = zeros(N, T);
    for k  = 1:N
        fitMFR(k,:) = exp([1 delt_fit(k,g) C_fit(k,1:THETA{g}(Lab(k)).p,g)]*...
            [THETA{g}(Lab(k)).mu ;ones(1,T) ;THETA{g}(Lab(k)).X]);
    end
    fitMFRTrace(:,:,g) = fitMFR;
    llhd_trace(g) = nansum(log(poisspdf(Y,fitMFR)), 'all')/...
        nansum(Y, 'all');
    
    
    % birth-death
    if g > bdstart
        [THETA{g}, C_fit(:,:,g)] = bdmcmc(THETA{g}, C_fit(:,:,g),...
            delt_fit(:,g), birth_rate,...
            birth_time,alpha, Y, unique(Lab), Lab, prior, p_max);
    end
    
    
    
    for kk = 1:nClus
        p_trace(kk,g) = THETA{g}(kk).p;
    end
    
    
    if g > bdstart
        figure(1)
        for kk = 1:nClus
            subplot(1,nClus,kk)
            histogram(p_trace(kk,bdstart:g))
            xline(2, 'r--', 'Linewidth', 2)
            title("p^{(g)} = "+p_trace(1,g))
        end
    end
    
    
    figure(2)
    hold on
    plot(delt_fit(:,g), 'k')
    plot(delta, 'r')
    hold off
    
    figure(3)
    for kk = 1:nClus
        subplot(1,nClus,kk)
        hold on
        plot(THETA{g}(kk).mu, 'k')
        plot(mu(kk,:), 'r', 'Linewidth', 2)
        hold off
        ylabel("racc = " + RACC(kk))
    end
    
    
    disp("iter " + g +": "+ toc)
end
toc









