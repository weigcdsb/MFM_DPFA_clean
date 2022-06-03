function x = myMvnrnd(mu, Sig)

% to debug
% mu = prior.BA0;
% Sig = theta.Q(k,k)*inv(prior.Lamb0);

y = randn(size(Sig,1),1);
R = chol(Sig,'lower');
x = mu + R*y;

end