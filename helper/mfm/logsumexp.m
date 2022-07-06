function  out = logsumexp(a,b)
m = max(a,b);
if m ==-Inf
    out = -Inf;
else
    out = log(exp(a-m) + exp(b-m)) + m;
end
end

