function lpdf = mvnlpdf(x, mu, invSig)

% -0.5*log(det(Sigma)) = 0.5*log(det(RR'))=
% 0.5*log(det(R'R)) = 0.5*log(det(R)^2) = log(det(R))
R = chol(invSig,'lower'); % sparse
dimx = length(x);
lpdf = -0.5*(x - mu)'*(R*R')*(x - mu) + sum(log(diag(R)))-0.5*dimx*log(2*pi);
end