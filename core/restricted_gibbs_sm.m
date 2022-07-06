function [zsb,tib,tjb,cib,cjb,log_p, ni, nj, iIdx, jIdx] =...
    restricted_gibbs_sm(zsa,zsb,tia,tib,tja,tjb,cia,cib,cja,cjb,...
    iIdx,jIdx, ni,nj,ism,jsm,S,ns,Y,b,prior, active, delt_s, C_s, R_all, useNormal)

% for debug
% zsa = zs;
% zsb = zs;
% tia = ti;
% tib = ti;
% tja = tj;
% tjb = tj;
% cia = ci;
% cib = ci;
% cja = cj;
% cjb = cj;
% active = true;

T = size(Y,2);
pia = tia.p;
pja = tja.p;

log_p = 0;
for ks = 1:ns
    k = S(ks);
    if k~= ism && k~=jsm
        if zsa(k) == cia
            ni = ni-1;
        else
            nj = nj-1;
        end
        
        
        lami = exp([1 delt_s(k) C_s(k,1:pia)]*...
            [tia.mu ;ones(1,T) ;tia.X]);
        Li = nansum(log(poisspdf(Y(k,:), lami)));
        lamj = exp([1 delt_s(k) C_s(k,1:pja)]*...
            [tja.mu ;ones(1,T) ;tja.X]);
        Lj = nansum(log(poisspdf(Y(k,:), lamj)));
        Pi = exp(log(ni+b)+Li - logsumexp(log(ni+b)+Li,log(nj+b)+Lj));
        
        
        if active
            if rand < Pi
                if zsa(k) == cja
                    jIdx = setdiff(jIdx, k);
                    iIdx = [iIdx; k];
                end
                zsb(k) = cib;
            else
                if zsa(k) == cia
                    iIdx = setdiff(iIdx, k);
                    jIdx = [jIdx; k];
                end
                zsb(k) = cjb;
            end
        end
        if zsb(k) == cib
            ni = ni+1;
            log_p = log_p + log(Pi);
        else
            nj = nj + 1;
            log_p = log_p + log(1-Pi);
        end
    end
end


prior.x0 = zeros(tia.p,1);
prior.Q0 = eye(tia.p);
prior.muC0 = zeros(tia.p,1);
prior.SigC0 = eye(tia.p);
C_s(iIdx,(tia.p+1):end) = 0;
[tib, log_pi] = update_factor_sm(tia,tib,delt_s(iIdx), C_s(iIdx,1:tia.p),Y(iIdx,:),...
    prior,active,true,R_all(iIdx,:),useNormal);

prior.x0 = zeros(tja.p,1);
prior.Q0 = eye(tja.p);
prior.muC0 = zeros(tja.p,1);
prior.SigC0 = eye(tja.p);
C_s(jIdx,(tja.p+1):end) = 0;
[tjb, log_pj] = update_factor_sm(tja,tjb,delt_s(jIdx), C_s(jIdx,1:tja.p),Y(jIdx,:),...
    prior,active,true,R_all(jIdx,:),useNormal);

log_p = log_p + log_pi + log_pj;


end