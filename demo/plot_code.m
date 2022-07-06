usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));

plotFolder_root = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot";

%% labeled
% load("C:\Users\gaw19004\Desktop\BDMCMC_PG\labeled_longBirthT.mat")
% plotFolder = plotFolder_root + "\labeled1";

load("C:\Users\gaw19004\Desktop\BDMCMC_PG\labeled_longBirthT_c2.mat")
plotFolder = plotFolder_root + "\labeled2";

% RACC: 0.5712
% RACC: 0.5327

idx = round(ng/4):ng;

% (1) trace of delta
delta_trace = figure;
plot(10:ng, delt_fit(:,10:ng)', '.')
title("trace for \delta_i")
xlabel("iteration")
set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
box off

set(delta_trace,'PaperUnits','inches','PaperPosition',[0 0 3 2])
saveas(delta_trace, plotFolder + "\1_delta.svg")
saveas(delta_trace, plotFolder + "\1_delta.png")


% (2) trace plot of ||mu||
muTrace = zeros(ng,1);
for gg = 1:ng
    muTrace(gg) = norm(THETA{gg}.mu,'fro');
end

mu_trace = figure;
plot(10:ng, muTrace(10:ng), '.')
title("trace for ||\mu||_2")
xlabel("iteration")
set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
box off

set(mu_trace,'PaperUnits','inches','PaperPosition',[0 0 3 2])
saveas(mu_trace, plotFolder + "\2_muNorm.svg")
saveas(mu_trace, plotFolder + "\2_muNorm.png")


% (3) fit of mu
muSum = zeros(nClus, T);

for gg = idx
    for kk = 1:nClus
        latidTmp = id2id(kk,p);
        muSum(kk,:) = muSum(kk,:) + THETA{gg}(kk).mu;
    end
end
mu_HPD = zeros(2,T,nClus);

for kk = 1:nClus
    mu_tmp = zeros(length(idx), T);
    
    c = 1;
    for gg = idx
        mu_tmp(c,:) = THETA{gg}(kk).mu;
        c = c+1;
    end
    mu_HPD(:,:,kk) = hpdi(mu_tmp, 95);
end

muMean = muSum/length(idx);

muFit = figure;
rrMax = 1;
cos_mu  = zeros(rrMax,1);
for rr = 1:rrMax
    subplot(rrMax,1,rr)
    hold on
    plot(mu(rr,:), 'k', 'LineWidth', 1)
    plot(muMean(rr,:), 'r', 'LineWidth', 1)
    plot(mu_HPD(1,:,rr), 'r--', 'LineWidth', 1)
    plot(mu_HPD(2,:,rr), 'r--', 'LineWidth', 1)
    xlabel('T')
    cos_mu(rr) = muMean(rr,:)*mu(rr,:)'/(norm(muMean(rr,:))*norm(mu(rr,:)));
    
    hold off
    set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
    box off
    if rr ==1
        title('\mu')
    end
end

% chain 1: cosine = 0.9641
% chain 2: cosine = 0.9686

set(muFit,'PaperUnits','inches','PaperPosition',[0 0 2 2])
saveas(muFit, plotFolder + "\3_mu.svg")
saveas(muFit, plotFolder + "\3_mu.png")


% (4) p
p_hist = figure;
histogram(p_trace(idx))
title("histogram of p")
set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
box off

set(p_hist,'PaperUnits','inches','PaperPosition',[0 0 2 2])
saveas(p_hist, plotFolder + "\4_p.svg")
saveas(p_hist, plotFolder + "\4_p.png")


%% unlabeled
clear all;close all;clc;
plotFolder_root = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot";

% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\unlabeled_chain2";
% plotFolder = plotFolder_root + "\unlabeled1";

% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\unlabeled_chain3";
% plotFolder = plotFolder_root + "\unlabeled2";

savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\unlabeled_chain4";
plotFolder = plotFolder_root + "\unlabeled3";

Z_all = [];

for kk = 1:10
    load(savedir + "\file" + kk + ".mat")
    Z_all = [Z_all Z_fit];
    clearvars -except Z_all t_trace savedir usr_dir plotFolder_root plotFolder Y Lab ng
end

% (1) trace of cluster number 
clusNum_trace = figure;
plot(2:ng, t_trace(2:end), '.')
xlabel('iteration')
title('# of clusters')
ylim([6 14])
set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
box off

set(clusNum_trace,'PaperUnits','inches','PaperPosition',[0 0 2 2])
saveas(clusNum_trace, plotFolder + "\1_trace.svg")
saveas(clusNum_trace, plotFolder + "\1_trace.png")

% (2) similarity matrix
idx = round(ng/4):ng;
nClus = 10;
n = 5;
N = size(Y,1);
simMat = zeros(N,N);
for g = idx
    for k = 1:size(simMat, 1)
        simMat(k,:) = simMat(k,:) + (Z_all(k, g) == Z_all(:, g))';
    end
end

sim_all = figure;
imagesc(simMat/length(idx))
hold on
for cc = 1:nClus
    yline(n*cc+0.5, 'k--', 'LineWidth', 1);
    xline(n*cc+0.5, 'k--', 'LineWidth', 1);
end
hold off
title('similarity matrix')
colormap(flipud(hot))
colorbar()
set(gca,'FontSize',10, 'LineWidth', 1.5,'TickDir','out')
box off

set(sim_all,'PaperUnits','inches','PaperPosition',[0 0 2.5 2])
saveas(sim_all, plotFolder + "\2_sim.svg")
saveas(sim_all, plotFolder + "\2_sim.png")

%% pixel data
clear all; close all;clc;

usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));
addpath(genpath(usr_dir + "\data"));

% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t1";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t1_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t2_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t3";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t3_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t4";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t4_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t5";
savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t5_v2";

Z_all = [];

for kk = 1:10
    load(savedir + "\file" + kk + ".mat")
    Z_all = [Z_all Z_fit];
    clearvars -except Z_all t_trace savedir usr_dir Y Lab ng clusIdx
end

save(savedir + "\plotData.mat");
% writematrix(Z_all, savedir + "\zLab.csv")


%
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t1";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t1_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t2_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t3";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t3_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t4";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t4_v2";
% savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t5";
savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t5_v2";


% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t1\v1";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t1\v2";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t2\v1";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t2\v2";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t3\v1";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t3\v2";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t4\v1";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t4\v2";
% plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t5\v1";
plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel\t5\v2";



load(savedir + "\plotData.mat")
zMaxPEAR = readmatrix(savedir + "\zMaxPEAR.csv");

N = size(Y,1);
zMaxPEAR = zMaxPEAR(:,2);
[zMAP_sort,id_all] = sort(zMaxPEAR);
id_sep = zeros(N,1);
for k = 1:length(unique(Lab))
    [~, idTmp] = sort(zMaxPEAR(Lab == k,:));
    id_sep(Lab == k) = sum(Lab < k) + idTmp;
end

idx = round(ng/4):ng;

% simMat
simMat = zeros(N,N);
for g = idx
    for k = 1:size(simMat, 1)
        simMat(k,:) = simMat(k,:) + (Z_all(k, g) == Z_all(:, g))';
    end
end

simPixel_all = figure;
imagesc(simMat(id_all,id_all)/length(idx))
colormap(flipud(hot))
colorbar()
set(gca,'FontSize',9, 'LineWidth', 1.5,'TickDir','out')
box off

set(simPixel_all,'PaperUnits','inches','PaperPosition',[0 0 2 1.5])
saveas(simPixel_all, plotFolder + "\1_simMat.svg")
saveas(simPixel_all, plotFolder + "\1_simMat.png")




% simMat: anatomy + maxPEAR
simPixel_sorted = figure;
imagesc(simMat(id_sep,id_sep)/length(idx))
hold on
nRegion = length(unique(Lab));
numEach = histcounts(Lab);

tickPos = zeros(1, nRegion);
for k = 1:nRegion
    yline(sum(numEach(1:k) + 0.5), 'b--', 'LineWidth', 2);
    xline(sum(numEach(1:k) + 0.5), 'b--', 'LineWidth', 2);
    tickPos(k) = sum(numEach(1:(k-1))) + numEach(k)/2;
end
yticks(tickPos)
yticklabels(clusIdx)
xticks(tickPos)
xticklabels(clusIdx)
hold off
colormap(flipud(hot))
colorbar()
set(gca,'FontSize',9, 'LineWidth', 1.5,'TickDir','out')
box off

set(simPixel_sorted,'PaperUnits','inches','PaperPosition',[0 0 2 1.5])
saveas(simPixel_sorted, plotFolder + "\2_simMat.svg")
saveas(simPixel_sorted, plotFolder + "\2_simMat.png")

%% ARI for pixel
clear all;close all;clc

plotFolder = "C:\Users\gaw19004\Documents\GitHub\MFM_DPFA_clean\plot\pixel";

ARI = readmatrix("C:\Users\gaw19004\Desktop\BDMCMC_PG\ari.csv");
ARI = ARI(:,2:end);

ari_heat = figure;
imagesc(tril(ARI))
colormap(flipud(gray(256)));
colorbar;
for r = 1:size(ARI,1)
    for c = 1:size(ARI,2)
        if(r>c)
            text(c,r,num2str(round(ARI(r,c),2)),'horizontalAlignment','center')
        elseif(r==c)
            text(c,r,num2str(round(ARI(r,c),2)),'horizontalAlignment','center','Color','w')
        end
    end
end 
set(gca,'FontSize',9, 'LineWidth', 1.5,'TickDir','out')
box off

cd(plotFolder)
set(ari_heat,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(ari_heat, 'ARI.svg')
saveas(ari_heat, 'ARI.png')

























