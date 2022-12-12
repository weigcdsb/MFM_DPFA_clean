usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));

savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\unlabeled_chain4";

%%
close all
rng(2)
n = 5;
nClus = 10;
N = n*nClus;
p = 2;
T = 1000;

Lab = repelem(1:nClus, n);

d = randn(n*nClus,1)*0.5;
C_trans = zeros(n*nClus, p*nClus);
for k = 1:nClus
    NtmpIdx = (n*(k-1)+1):(n*k);
    d(NtmpIdx) = d(NtmpIdx) - mean(d(NtmpIdx));
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
        %         disp(k)
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
logLam(:,1) = d + C_trans*X(:,1) + repelem(mu(:,1), n);


for t=2:T
    logLam(:, t) = d + C_trans*X(:,t) + repelem(mu(:,t), n);
end

%%
Y = poissrnd(exp(logLam));
figure(1)
subplot(1,2,1)
imagesc(exp(logLam))
colorbar()
subplot(1,2,2)
imagesc(exp(mu))
colorbar()

%%
rng(4)
close all
ng = 10000;
burnIn = round(ng/10);
nSepFile = 10;
iterEach = ng/nSepFile;

p_max = 20;

t_max = N;
lAbsGam = @(x) log(abs(gamma(x)));

% DPMM:
% DPMM = true;
% alpha_random = true;
sigma_alpha = 0.1; % scale for MH proposals in alpha move
alphaDP = 1;
% log_v = (1:t_max+1)*log(alphaDP) - lAbsGam(alphaDP+N) + lAbsGam(alphaDP);
% a = 1;
% b = 0;

% MFM:
DPMM = false;
alpha_random = false;
MFMgamma = 1;
% K ~ Geometric(r)
r = 0.2;
log_pk = @(k) log(r) + (k-1)*log(1-r);
a = MFMgamma;
b = MFMgamma;
log_v = MFMcoeff(log_pk, MFMgamma, N, t_max + 1);

logNb = log((1:N) + b);

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
t_fit = zeros(1, iterEach);
Z_fit = zeros(N, iterEach);
numClus_fit = zeros(t_max + 3, iterEach);

% start from 1 cluster
t_fit(1) = 1;
Z_fit(:,1) = ones(N, 1);
numClus_fit(1,1) = N;
actList = zeros(t_max+3,iterEach); actList(1,1) = 1;
c_next = 2;

delt_fit = zeros(N,iterEach);
C_fit = zeros(N,p_max,iterEach);

delt_fit(:,1) = normrnd(prior.delt0, prior.Sigdelt0, N,1);

for k = 1:t_fit(1)
    
    c = actList(k,1);
    idx_tmp = find(Z_fit(:,1) == c);
    N_tmp = length(idx_tmp);
    p_tmp = 1;
    
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
fitMFRTrace = zeros(N, T, iterEach);
llhd_trace = zeros(iterEach,1);

% tuning parameter: NB dispersion
R_all = 10*ones(N,T);

% BDMCMC
birth_rate = .5;
birth_time = 1;
alpha = 2;

%% MCMC
if isempty(gcp('nocreate'))
    poolobj = parpool('IdleTimeout',120);
end


useSplitMerge = true;
splitPeriod = 10;

simMat = zeros(N,N);
for k = 1:size(simMat, 1)
    simMat(k,:) = simMat(k,:) + (Z_fit(k, 1) == Z_fit(:, 1))';
end

% AGS
bdstart = 2;% burnIn/2;
b0 = 0.1;
b1 = 1e-4;

simMat = zeros(N,N);
for k = 1:size(simMat, 1)
    simMat(k,:) = simMat(k,:) + (Z_fit(k, 1) == Z_fit(:, 1))';
end

p_trace = ones(N, iterEach);
nRep = 4;


t_trace = zeros(1, ng);
t_trace(1) = t_fit(1);

ct = 1;
tic;
for ff = 1:nSepFile
    for g = 2:iterEach
        ct = ct + 1;
        
        THETA{g} = THETA{g-1};
        delt_fit(:,g) = delt_fit(:,g-1);
        C_fit(:,:,g) = C_fit(:,:,g-1);
        llhd_pre = llhd_trace(g-1);
        
        Z_fit(:,g) = Z_fit(:,g-1);
        numClus_fit(:,g) = numClus_fit(:,g-1);
        t_fit(g) = t_fit(g-1);
        actList(:,g) = actList(:,g-1);
        
        for rr = 1:nRep
            
            % (1) update cluster parameters (PG+FFBS+MH)
            [THETA{g}, ~] = update_clusParam_all(THETA{g}, C_fit(:,:,g),...
                delt_fit(:,g), unique(Z_fit(:,g)),...
                Z_fit(:,g), Y, prior, R_all, false);
            
            % (2) update non-cluster parameters: intercept + loading (HMC)
            [delt_fit(:,g), C_fit(:,:,g)] =...
                update_deltC_HMC(delt_fit(:,g), C_fit(:,:,g),...
                THETA{g},Z_fit(:,g),Y,prior);
            
            % (3) projection
            [THETA{g},  C_fit(:,:,g), delt_fit(:,g)] = constraintProj(THETA{g},...
                C_fit(:,:,g), delt_fit(:,g),unique(Z_fit(:,g)), Z_fit(:,g));
            
        end
        
        fitMFR = zeros(N, T);
        for k  = 1:N
            fitMFR(k,:) = exp([1 delt_fit(k,g) C_fit(k,1:THETA{g}(Z_fit(k,g)).p,g)]*...
                [THETA{g}(Z_fit(k,g)).mu ;ones(1,T) ;THETA{g}(Z_fit(k,g)).X]);
        end
        fitMFRTrace(:,:,g) = fitMFR;
        llhd_trace(g) = nansum(log(poisspdf(Y,fitMFR)), 'all')/...
            nansum(Y, 'all');
        
        % (4) update cluster...
        
        % use split-merge?
        if useSplitMerge && (mod(g,splitPeriod) == 0)
            [Z_fit(:,g), numClus_fit(:,g), t_fit(g), actList(:,g),...
                THETA{g},delt_fit(:,g),C_fit(:,:,g)] =...
                splitMerge(Z_fit(:,g), numClus_fit(:,g), t_fit(g), actList(:,g),...
                THETA{g},delt_fit(:,g),C_fit(:,:,g),Y,R_all, prior,...
                a, b, log_v,binornd(1,.5), 5, 5); % binornd(1,.5); nan;
            c_next = ordered_next(actList(:,g));
        end
        
        [Z_fit(:,g), numClus_fit(:,g), t_fit(g), actList(:,g), c_next, THETA{g}] =...
            update_cluster(Z_fit(:,g), numClus_fit(:,g), t_fit(g), actList(:,g), c_next,...
            DPMM, alpha_random, a, log_v, logNb,...
            THETA{g}, delt_fit(:,g), Y, prior,alphaDP,sigma_alpha,t_max);
        
        for j = 1:t_fit(g)
            c = actList(j,g);
            obsIdx = find(Z_fit(:,g) == c);
            C_add = normrnd(0,1, length(obsIdx), THETA{g}(c).p);
            C_ori = C_fit(obsIdx, 1:THETA{g}(c).p,g);
            C_ori(C_ori == 0) = C_add(C_ori == 0);
            
            C_fit(obsIdx, 1:THETA{g}(c).p,g) = C_ori;
        end
        
        % (5) update p: BDMCMC
        % birth-death
        if g > bdstart % && binornd(1, exp(-b0 - b1*g)) == 1
            [THETA{g}, C_fit(:,:,g)] = bdmcmc(THETA{g}, C_fit(:,:,g),...
                delt_fit(:,g), birth_rate,...
                birth_time,alpha, Y, unique(Z_fit(:,g)), Z_fit(:,g), prior, p_max);
        end
        
        
        t_trace(ct) = t_fit(g);
        disp("file = " + ff + ", iter = " + g + ": " + toc)
        
        figure(1)
        clusterPlot(Y, Z_fit(:,g)')
        
        figure(2)
        plot(t_trace(1:ct))
        
        
    end
    
    save(savedir+"\file"+ff+".mat");
    
    % one more step to store in 1st iteration
    ct = ct + 1;
    
    THETA{1} = THETA{iterEach};
    delt_fit(:,1) = delt_fit(:,iterEach);
    C_fit(:,:,1) = C_fit(:,:,iterEach);
    llhd_pre = llhd_trace(iterEach);
    
    Z_fit(:,1) = Z_fit(:,iterEach);
    numClus_fit(:,1) = numClus_fit(:,iterEach);
    t_fit(1) = t_fit(iterEach);
    actList(:,1) = actList(:,iterEach);
    
    
    for rr = 1:nRep
        
        % (1) update cluster parameters (PG+FFBS+MH)
        [THETA{1}, ~] = update_clusParam_all(THETA{1}, C_fit(:,:,1),...
            delt_fit(:,1), unique(Z_fit(:,1)),...
            Z_fit(:,1), Y, prior, R_all, false);
        
        % (2) update non-cluster parameters: intercept + loading (HMC)
        [delt_fit(:,1), C_fit(:,:,1)] =...
            update_deltC_HMC(delt_fit(:,1), C_fit(:,:,1),...
            THETA{1},Z_fit(:,1),Y,prior);
        
        % (3) projection
        [THETA{1},  C_fit(:,:,1), delt_fit(:,1)] = constraintProj(THETA{1},...
            C_fit(:,:,1), delt_fit(:,1),unique(Z_fit(:,1)), Z_fit(:,1));
        
    end
    
    fitMFR = zeros(N, T);
    for k  = 1:N
        fitMFR(k,:) = exp([1 delt_fit(k,1) C_fit(k,1:THETA{1}(Z_fit(k,1)).p,1)]*...
            [THETA{1}(Z_fit(k,1)).mu ;ones(1,T) ;THETA{1}(Z_fit(k,1)).X]);
    end
    fitMFRTrace(:,:,1) = fitMFR;
    llhd_trace(1) = nansum(log(poisspdf(Y,fitMFR)), 'all')/...
        nansum(Y, 'all');
    
    % (4) update cluster...
    % use split-merge?
    if useSplitMerge && (mod(1+ff*iterEach,splitPeriod) == 0)
        [Z_fit(:,1), numClus_fit(:,1), t_fit(1), actList(:,1),...
            THETA{1},delt_fit(:,1),C_fit(:,:,1)] =...
            splitMerge(Z_fit(:,1), numClus_fit(:,1), t_fit(1), actList(:,1),...
            THETA{1},delt_fit(:,1),C_fit(:,:,1),Y,R_all, prior,...
            a, b, log_v,binornd(1,.5), 5, 5); % binornd(1,.5); nan;
        c_next = ordered_next(actList(:,1));
    end
    
    [Z_fit(:,1), numClus_fit(:,1), t_fit(1), actList(:,1), c_next, THETA{1}] =...
        update_cluster(Z_fit(:,1), numClus_fit(:,1), t_fit(1), actList(:,1), c_next,...
        DPMM, alpha_random, a, log_v, logNb,...
        THETA{1}, delt_fit(:,1), Y, prior,alphaDP,sigma_alpha,t_max);
    
    for j = 1:t_fit(1)
        c = actList(j,1);
        obsIdx = find(Z_fit(:,1) == c);
        C_add = normrnd(0,1, length(obsIdx), THETA{1}(c).p);
        C_ori = C_fit(obsIdx, 1:THETA{1}(c).p,1);
        C_ori(C_ori == 0) = C_add(C_ori == 0);
        
        C_fit(obsIdx, 1:THETA{1}(c).p,1) = C_ori;
    end
    
    % (5) update p: BDMCMC
    % birth-death
    if 1+ff*iterEach > bdstart % && binornd(1, exp(-b0 - b1*(1+ff*iterEach))) == 1
        [THETA{1}, C_fit(:,:,1)] = bdmcmc(THETA{1}, C_fit(:,:,1),...
            delt_fit(:,1), birth_rate,...
            birth_time,alpha, Y, unique(Z_fit(:,1)), Z_fit(:,1), prior, p_max);
    end
    
    t_trace(ct) = t_fit(1);
    disp("file = " + (ff + 1) + ", iter = 1: " + toc)
    

    figure(1)
    clusterPlot(Y, Z_fit(:,1)')
    
    figure(2)
    plot(t_trace(1:ct))
    
    actTmp = unique(actList(:,1));
    if sum(actTmp(actTmp~=0) ~= sort(unique(Z_fit(:,1)))) ~= 0
        error('mismatch')
    end
    
end


% check_code

% tab = tabulate(Z_fit(:,g));
% test = zeros(t_max + 3,1);
% test(tab(:,1)) = tab(:,2);
% if sum(numClus_fit(:,g) ~= test) > 0
%     eflg = 1;
%     error('mismatch')
% end