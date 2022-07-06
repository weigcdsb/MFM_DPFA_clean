function j = randlogp(log_p, k)
log_s = -Inf;
for j = 1:k
    log_s = logsumexp(log_s,log_p(j)); 
end
p = log_p;
for j = 1:k
    p(j) = exp(log_p(j)-log_s); 
end
j = randp(p,k);
end

function j = randp(p,k)
s = 0;
for j = 1:k
    s = s + p(j);
end
u = rand()*s;
j = 1;
C = p(1);
while u > C
    j = j + 1;
    C = C + p(j);
end
end










