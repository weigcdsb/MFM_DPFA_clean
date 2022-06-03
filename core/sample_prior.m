function theta =...
    sample_prior(prior, T, p, sampleDynamics)


% linear dynamics: b, A & Q
theta.A = eye(p+1);
theta.b = zeros(p+1,1);
theta.Q = eye(p+1)*prior.sig20;

% mean & latent: d, X
muX = ones(1+p,T)*Inf;
muX(:,1) = myMvnrnd([prior.mu0;prior.x0], blkdiag(prior.Sigmu0, prior.Q0));


if sampleDynamics
    for k = 1:p+1
        theta.Q(k,k) = 1/gamrnd(prior.nu0/2,...
            2/(prior.nu0*prior.sig20));
        baSamp = myMvnrnd(prior.BA0,...
            theta.Q(k,k)*inv(prior.Lamb0));
        
        theta.b(k) = baSamp(1);
        theta.A(k,k) = baSamp(2);
    end
end

for i=1:(p+1)
    while(sum(muX(i,:) > 2) > 1)
        k = ceil(rand()*25)+10;
        muX(i,:) = interp1(linspace(0,1,k),randn(k,1)*0.3,linspace(0,1,T),'spline');
    end
end


theta.mu = muX(1,:);
theta.X = muX(2:end,:);

muBar = mean(theta.mu);
theta.mu = theta.mu - muBar;

XBar = mean(theta.X, 2);
theta.X = theta.X - XBar;


end