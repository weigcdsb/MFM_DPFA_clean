usr_dir = "C:\Users\gaw19004\Documents\GitHub";
addpath(genpath(usr_dir + "\MFM_DPFA_clean"));
addpath(genpath(usr_dir + "\data"));

savedir = "C:\Users\gaw19004\Desktop\BDMCMC_PG\pixel_t1";

%% load data
Z_all = [];

for kk = 1:10
    load(savedir + "\file" + kk + ".mat")
    Z_all = [Z_all Z_fit];
    clearvars -except Z_all t_trace savedir usr_dir Y Lab ng clusIdx
end


%%
% trace of cluster number
plot(2:ng, t_trace(2:ng), '.')

% similarity matrix

% raw
N = size(Y,1);
simMat = zeros(N,N);
for g = 1:ng
    for k = 1:size(simMat, 1)
        simMat(k,:) = simMat(k,:) + (Z_all(k, g) == Z_all(:, g))';
    end
end

imagesc(simMat/ng)
colormap(flipud(hot))
colorbar()

% writematrix(Z_all, savedir + "\zLab.csv")
zMAP = readmatrix(savedir + "\zMaxPEAR.csv");
zMAP = zMAP(:,2);
[zMAP_sort,id_all] = sort(zMAP);
id_sep = zeros(N,1);
for k = 1:length(unique(Lab))
    [~, idTmp] = sort(zMAP(Lab == k,:));
    id_sep(Lab == k) = sum(Lab < k) + idTmp;
end

idx = round(ng/4):ng;

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





