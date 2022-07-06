function log_v = MFMcoeff(log_pk, MFMgamma, n, upto)

lAbsGam = @(x) log(abs(gamma(x)));
tol = 1e-12;
log_v = zeros(upto,1);

for t = 1:upto
    if t > n
        log_v(t) = -Inf;
        continue;
    end
    [a,c,k,p] = deal(0, -Inf,1,0);
    while(abs(a-c) > tol || p < 1- tol)
        if k>= t
            a=c;
            b = lAbsGam(k+1)-lAbsGam(k-t+1) -...
                lAbsGam(k*MFMgamma+n)+lAbsGam(k*MFMgamma) + log_pk(k);
            c = logsumexp(a,b);
        end
        p = p + exp(log_pk(k));
        k = k+1;
    end
    log_v(t) = c;
end

end