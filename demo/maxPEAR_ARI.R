library(label.switching)
library(mclust)
library(mcclust)

folder_root <- "C:/Users/gaw19004/Desktop/BDMCMC_PG"


for(i in 1:5){
  for(j in 1:2){
    if(j == 1){
      folder_tmp <- paste0(folder_root, "/pixel_t", i)
    }else{
      folder_tmp <- paste0(folder_root, "/pixel_t", i, "_v2")
    }
    
    z <- read.csv(paste0(folder_tmp, "/zLab.csv"),header = F)
    
    ng <- dim(z)[2]
    idx <- round(ng/4):ng
    z_use <- t(as.matrix(z[,idx]))
    
    psm <- comp.psm(z_use)
    res <- maxpear(psm)
    z2 <- res$cl
    write.csv(z2,paste0(folder_tmp, "/zMaxPEAR.csv"))
    
    
  }
}

######################################################
#### calculate ARI

zlab_each <- matrix(NA, 5, length(z2))
zlab <- list(zlab_each, zlab_each)

for(ep in 1:5){
  for(rep in 1:2){
    if(rep == 1){
      folder_tmp <- paste0(folder_root, "/pixel_t", ep)
    }else{
      folder_tmp <- paste0(folder_root, "/pixel_t", ep, "_v2")
    }
    zlab[[rep]][ep,] <- read.csv(paste0(folder_tmp, "/zMaxPEAR.csv"))[,-1]
  }
}

ari_self <- rep(NA, 5)
for(ep in 1:5){
  ari_self[ep] <- adjustedRandIndex(zlab[[1]][ep,], zlab[[2]][ep,])
}

ari_cross <- array(NA, dim = c(5,5,4))
count <- 1
for(rep1 in 1:2){
  for(rep2 in 1:2){
    for(r in 2:5){
      for(c in 1:(r-1)){
        ari_cross[r,c,count] <- adjustedRandIndex(zlab[[rep1]][r,], zlab[[rep2]][c,])
      }
    }
    count <- count + 1
  }
}


ARImean <- apply(ari_cross, c(1,2), mean)
diag(ARImean) <- ari_self
ARImin <- apply(ari_cross, c(1,2), min)
ARImax <- apply(ari_cross, c(1,2), max)

write.csv(ARImean,"C:/Users/gaw19004/Desktop/BDMCMC_PG/ari.csv")
write.csv(ARImin,"C:/Users/gaw19004/Desktop/BDMCMC_PG/ariMin.csv")
write.csv(ARImax,"C:/Users/gaw19004/Desktop/BDMCMC_PG/ariMax.csv")





