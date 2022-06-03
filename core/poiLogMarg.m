function logMar = poiLogMarg(Yreg, Xreg, offset)

% to debug
% Yreg = Y(ii,:)';
% Xreg = THETA{g}(c_prop).X';
% offset = delta_tmp(ii)*ones(T,1) + THETA_b(cc).mu';

% way 1: do gamma approximation to lambda
T = length(Yreg);
p = size(Xreg, 2);
% XX = diag(Xreg*Xreg');
% for kk = 1:T
%     f = @(lam) poisspdf(Yreg(kk), lam).*lognpdf(lam,offset(kk),XX(kk));
%     logMarvec(kk) = log(integral(f,0,Inf));
% end

if p ~=0
    XX = diag(Xreg*Xreg');
    avec = XX.^(-1);
    bvec = XX.*exp(offset);
    probvec = 1./(1 + bvec);
    logMarvec = log(nbinpdf(Yreg, avec, probvec));
    logMar =nansum(logMarvec);
else
    logMar =nansum(-exp(offset) + Yreg.*offset);
end




% way 2: Laplace approximation...




end