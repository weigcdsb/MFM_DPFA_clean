function [Z_b, numClus_b, t_b, actList, c_next, THETA_b] =...
    update_cluster(Z_a, numClus_a, t_a, actList, c_next,...
    DPMM, alpha_random, a, log_v, logNb,...
    THETA_a, delta_tmp, Y_tmp, prior, alphaDP,sigma_alpha,t_max)

% to debug
% Z_a = Z_fit(:,g);
% numClus_a = numClus_fit(:,g);
% t_a = t_fit(g);
% THETA_a = THETA{g};
% delta_tmp = delt_fit(:,g);

N = size(Y_tmp,1);
T = size(Y_tmp,2);

lAbsGam = @(x) log(abs(gamma(x)));
Z_b = Z_a;
numClus_b = numClus_a;
t_b = t_a;
THETA_b = THETA_a;

if DPMM && alpha_random
    % MH move for DP concentration parameter (using p_alpha(a) = exp(-a) = Exp(a|1))
    aprop = alphaDP*exp(rand*sigma_alpha);
    top = t_a*log(aprop) - lAbsGam(aprop+N) + lAbsGam(aprop) - aprop + log(aprop);
    bot = t_a*log(alphaDP) - lAbsGam(alphaDP+N) +...
        lAbsGam(alphaDP) - alphaDP + log(alphaDP);
    if rand < min(1, exp(top-bot))
        alphaDP = aprop;
    end
    log_v = (1:t_max+1)*log(alphaDP) - lAbsGam(alphaDP+N) + lAbsGam(alphaDP);
end

for ii = 1:N
    
    % (a) remove point i from its cluster
    c = Z_b(ii);
    numClus_b(c) = numClus_b(c) - 1;
    if(numClus_b(c) > 0)
        c_prop = c_next;
        prior.x0 = 0;
        prior.Q0 = 1;
        prior.muC0 = 0;
        prior.SigC0 = 1;
        theta_tmp = sample_prior(prior, T, 1, false);
        theta_tmp.p = 1;
        THETA_b(c_prop) = theta_tmp;
    else
        c_prop = c;
        actList = ordered_remove(c, actList, t_b);
        t_b = t_b - 1;
    end
    
    %(b) compute probabilities for resampling
    log_p = zeros(t_b+1,1);
    for j = 1:t_b
        cc = actList(j);
        logMar = poiLogMarg(Y_tmp(ii,:)', THETA_b(cc).X',...
            delta_tmp(ii)*ones(T,1) + THETA_b(cc).mu');
        
        log_p(j) = logNb(numClus_b(cc)) + logMar;
    end
    
    logMar = poiLogMarg(Y_tmp(ii,:)', THETA_b(c_prop).X',...
        delta_tmp(ii)*ones(T,1) + THETA_b(c_prop).mu');
    
    
    log_p(t_b+1) = log_v(t_b+1)-log_v(t_b) +...
        log(a) + logMar;
    
    % (c) sample a new cluster for it
    j = randlogp(log_p, t_b+1);
    
    % (d) add point i to its new clusters
    if j <= t_b
        c = actList(j);
    else
        c = c_prop;
        actList = ordered_insert(c, actList, t_b);
        t_b = t_b + 1;
        c_next = ordered_next(actList);
    end
    
    Z_b(ii) = c;
    numClus_b(c) = numClus_b(c) + 1;
end


end