function latID = id2id(clusID, p)

nClus_tmp = length(clusID);
latID = zeros(nClus_tmp*p, 1);
for k = 1:nClus_tmp
    cid = clusID(k);
    latID(((k-1)*p+1):(k*p)) = ((cid-1)*p + 1):(cid*p);
end

end