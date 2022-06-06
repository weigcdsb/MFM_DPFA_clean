function [delta_b, C_b] = update_deltC_HMC(delta_a, C_a, THETA_a,Z_tmp, Y_tmp, prior)

% for debugging
% delta_a = delt_fit(:,g-1);
% C_a = C_fit(:,:,g-1);
% THETA_a = THETA{g};
% Z_tmp = Lab;
% Y_tmp = Y;

N = size(Y_tmp,1);
T = size(Y_tmp,2);

delta_b = delta_a;

C_b_cell = cell(N,1);
C_b_cell(:) = {zeros(1,size(C_a,2))};

parfor ii = 1:N
    
    idxTmp =  find(~isnan(Y_tmp(ii,:)));
    X_tmp = [ones(length(idxTmp),1) THETA_a(Z_tmp(ii)).X(:,idxTmp)'];
    lamdeltC = @(deltC) exp(THETA_a(Z_tmp(ii)).mu(idxTmp)' + X_tmp*deltC);
    
    prior_tmp = prior;
    
    prior_tmp.muC0 = zeros(THETA_a(Z_tmp(ii)).p,1);
    prior_tmp.SigC0 = eye(THETA_a(Z_tmp(ii)).p);
    
    lpdf = @(deltC) sum(-lamdeltC(deltC) +...
        Y_tmp(ii,idxTmp)'.*log(lamdeltC(deltC) + (lamdeltC(deltC) == 0))) +...
        log(mvnpdf(deltC, [prior_tmp.delt0;prior_tmp.muC0],...
        blkdiag(prior_tmp.Sigdelt0, prior_tmp.SigC0)));
    glpdf = @(deltC) X_tmp'*(Y_tmp(ii,idxTmp)' - lamdeltC(deltC)) -...
        blkdiag(prior_tmp.Sigdelt0, prior_tmp.SigC0)\(deltC - [prior_tmp.delt0;prior_tmp.muC0]);
    fg=@(deltC)logpdf(deltC, lpdf, glpdf);
    
    deltC0 = [delta_a(ii); C_a(ii,1:THETA_a(Z_tmp(ii)).p)'];
    smp = hmcSampler(fg,deltC0);
    samp = drawSamples(smp,'Burnin',0,'NumSamples',1,'StartPoint',deltC0);
    
    delta_b(ii) = samp(1);
    C_b_cell{ii}(1:THETA_a(Z_tmp(ii)).p) = samp(2:end);
end

C_b = cell2mat(C_b_cell);

end