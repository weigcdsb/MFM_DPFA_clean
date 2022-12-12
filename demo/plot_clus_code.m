usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));
savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\unlabeled_chain4";
plotFolder_root = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\unlabeled_clus_param";

%% load results
ng = 10000;
T = 1000;
z_clus_plot = zeros(3,ng);
p_clus_plot = zeros(3,ng);
for qq = 1:3
    mu_clus_plot{qq} = zeros(ng, T);
end

for kk = 1:10
    load(savedir + "\file" + kk + ".mat")
    
    ng_chunk = size(Z_fit, 2);
    T = size(Y,2);
    iter_idx = ((kk-1)*ng_chunk+1):(kk*ng_chunk);
    
    % cluster 1-3
    for qq = 1:3
        z_clus_plot(qq,iter_idx) = mode(Z_fit(((qq-1)*5+1):(qq*5),:), 1);
    end
    
    for gg = 1:ng_chunk
        for qq = 1:3
            mu_clus_plot{qq}(iter_idx(gg),:) = THETA{gg}(z_clus_plot(qq,iter_idx(gg))).mu;
            p_clus_plot(qq,iter_idx(gg)) = THETA{gg}(z_clus_plot(qq,iter_idx(gg))).p;
        end
    end
    
end

%% plotting...
idx = round(ng/4):ng;

% (1) histogram of p

cd(plotFolder_root + "\p")
for qq = 1:3
    trace_p = figure;
    plot(p_clus_plot(qq,:))
    yline(2, 'r--', 'LineWidth', 2)
    title("trace for p")
    xlabel("iteration")
    set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(trace_p,'PaperUnits','inches','PaperPosition',[0 0 2 2])
    saveas(trace_p, qq + "_trace_p.svg")
    saveas(trace_p, qq + "_trace_p.png")
    
    hist_p = figure;
    histogram(p_clus_plot(qq,idx))
    title("histogram for p")
    xline(2, 'r--', 'LineWidth', 2)
    xlabel("p")
    xlim([0 max(p_clus_plot, [], 'all')])
    set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(hist_p,'PaperUnits','inches','PaperPosition',[0 0 2 2])
    saveas(hist_p, qq + "_hist_p.svg")
    saveas(hist_p, qq + "_hist_p.png")
    
end


cd(plotFolder_root + "\mu")
% (2) trace of mu norm
muNorm = zeros(3, ng);
for gg = 1:ng
    for qq = 1:3
        muNorm(qq,gg) = norm(mu_clus_plot{qq}(gg,:), 2);
    end
end

for qq = 1:3
    mu_trace = figure;
    plot(10:ng, muNorm(qq,10:ng))
    yline(norm(mu(qq,:), 2), 'r--', 'LineWidth', 2)
    title("trace for ||\mu||_2")
    xlabel("iteration")
    set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(mu_trace,'PaperUnits','inches','PaperPosition',[0 0 2 2])
    saveas(mu_trace, qq + "_muNorm.svg")
    saveas(mu_trace, qq + "_muNorm.png")
    
end

% (3) fit of mu
mu_HPD = zeros(2,T,3);
muMean = zeros(3, T);
for qq = 1:3
    mu_HPD(:,:,qq) = hpdi(mu_clus_plot{qq}(idx,:), 95);
    muMean(qq,:) = mean(mu_clus_plot{qq}(idx,:), 1); 
end

cos_mu  = zeros(3,1);
for qq = 1:3
    
    muFit = figure;
    cos_mu(qq) = muMean(qq,:)*mu(qq,:)'/(norm(muMean(qq,:))*norm(mu(qq,:)));
    hold on
    plot(mu(qq,:), 'k', 'LineWidth', 1)
    plot(muMean(qq,:), 'r', 'LineWidth', 1)
    plot(mu_HPD(1,:,qq), 'r--', 'LineWidth', 1)
    plot(mu_HPD(2,:,qq), 'r--', 'LineWidth', 1)
    hold off
    title('fit of \mu')
    xlabel('T')
    set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(muFit,'PaperUnits','inches','PaperPosition',[0 0 2 2])
    saveas(muFit, qq + "_muFit.svg")
    saveas(muFit, qq + "_muFit.png")
    
end



















