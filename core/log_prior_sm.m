function lpdf = log_prior_sm(theta, prior, delt, C)

% to debug
% theta = ti;
% delt = delt_s(mIdx,:);
% C = C_s(mIdx,:);

% mudc: [3×1 double]
% Sigdc: [3×3 double]
% d: [30×1 double]
% C: [30×2 double]
% A: [2×2 double]
% b: [2×1 double]
% Q: [2×2 double]
% Xori: [2×1000 double]
% x0: [2×1 double]
% X: [2×1000 double]

p = theta.p;
T = size(theta.X, 2);

prior.x0 = zeros(p,1);
prior.Q0 = eye(p);
prior.muC0 = zeros(p,1);
prior.SigC0 = eye(p);
lpdf = 0;

% (1) muX
muX =  [theta.mu;theta.X];
lpdf = lpdf + mvnlpdf(muX(:,1), [prior.mu0;prior.x0], inv(blkdiag(prior.Sigmu0, prior.Q0)));

% for t = 2:T
%     lpdf = lpdf + mvnlpdf(muX(:,t), theta.b + theta.A*muX(:, t-1), inv(theta.Q));
% end


% (2) Q, b & A
for kk = 1:(p+1)
    lpdf = lpdf + log(gampdf(1/theta.Q(kk,kk), prior.nu0/2,...
        2/(prior.nu0*prior.sig20)));
    baTmp = [theta.b(kk) theta.A(kk,kk)]';
    lpdf = lpdf + mvnlpdf(baTmp, prior.BA0, theta.Q(kk,kk)*inv(prior.Lamb0));
end


% (3) delta & C
prior.muC0 = zeros(p,1);
prior.SigC0 = eye(p);

deltC = [delt C(:, 1:p)];
for jj = 1:size(deltC,1)
   lpdf = lpdf + ...
   log(mvnpdf(deltC(jj,:)', [prior.delt0;prior.muC0],...
        blkdiag(prior.Sigdelt0, prior.SigC0))); 
end


end