function logMar = poiLogMarg(Yreg, Xreg, offset)

% to debug
% Yreg = Y(ii,:)';
% Xreg = THETA_b(cc).X';
% offset = delta_tmp(ii)*ones(T,1) + THETA_b(cc).mu';

% Yreg = reshape(Y(obsIdx,:)', [], 1);
% Xreg = repmat(THETA{g}(c).X', N_tmp,1);
% offset = kron(delt_fit(obsIdx,g), ones(T,1)) + repmat(THETA{g}(c).mu', N_tmp,1);

% way 1: do gamma approximation to lambda
T = length(Yreg);
p = size(Xreg, 2);
% XX = diag(Xreg*Xreg');
% for kk = 1:T
%     f = @(lam) poisspdf(Yreg(kk), lam).*lognpdf(lam,offset(kk),XX(kk));
%     logMarvec(kk) = log(integral(f,0,Inf));
% end



if p ~=0
    XX = sum(Xreg.^2,2);
    avec = XX.^(-1);
    bvec = XX.*exp(offset);
    probvec = 1./(1 + bvec);
    logMarvec = log(nbinpdf(Yreg, avec, probvec));
    logMar =nansum(logMarvec);
else
    logMarvec = log(poisspdf(Yreg, exp(offset)));
    logMar =nansum(logMarvec);
end




% way 2: Laplace approximation...




end